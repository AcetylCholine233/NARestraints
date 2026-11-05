from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio.PDB import PDBParser

from .builder import first_model, pair_residue_from_pdb
from .phenix import (
    AT_ANGLE_RESTRAINTS,
    AT_BONDS,
    AT_PARALLEL_A,
    AT_PARALLEL_T,
    GC_ANGLE_RESTRAINTS,
    GC_BONDS,
    GC_PARALLEL_C,
    GC_PARALLEL_G,
    DT_ANGLE_RESTRAINTS,
    DT_BONDS,
    PairResidue,
    atom_name,
)
from .metal_pairs import (
    DEFAULT_MAX_DISTANCE,
    MetalMediatedPair,
    find_metal_mediated_pairs,
    metal_bound_residues,
)
from .recipe_library import recipe_for
from .validator import validate_and_fix_pdb
from .residue_library import load_residue_records


# Deliberately forgiving because this is intended for low-resolution models.
SEARCH_RADIUS = 8.0
BOND_SOFT_SIGMA = 0.40
MAX_PLANE_ANGLE = 25.0
MIN_SCORE = 45.0
MIN_RUN_SCORE = 35.0


@dataclass(frozen=True)
class Candidate:
    first: PairResidue
    second: PairResidue
    recipe: str
    score: float
    bond_distances: tuple[float, ...]
    angle_errors: tuple[float, ...]
    plane_angle: float
    noncanonical: bool = False


def _template(recipe_name: str):
    if recipe_name == "AT":
        return AT_ANGLE_RESTRAINTS, AT_BONDS, AT_PARALLEL_A, AT_PARALLEL_T
    if recipe_name == "D_T":
        return DT_ANGLE_RESTRAINTS, DT_BONDS, GC_PARALLEL_G, GC_PARALLEL_C
    if recipe_name == "AG_IX":
        from .phenix import AG_IX_ANGLES, AG_IX_BONDS
        return AG_IX_ANGLES, AG_IX_BONDS, GC_PARALLEL_G, GC_PARALLEL_G
    if recipe_name == "GU_XXVIII":
        from .phenix import GU_XXVIII_ANGLES, GU_XXVIII_BONDS
        return GU_XXVIII_ANGLES, GU_XXVIII_BONDS, GC_PARALLEL_G, AT_PARALLEL_T
    return GC_ANGLE_RESTRAINTS, GC_BONDS, GC_PARALLEL_G, GC_PARALLEL_C


def _role_residues(first: PairResidue, second: PairResidue, recipe):
    if (
        first.base_class in recipe.role1_categories
        and second.base_class in recipe.role2_categories
    ):
        return first, second
    if (
        second.base_class in recipe.role1_categories
        and first.base_class in recipe.role2_categories
    ):
        return second, first
    return None


def _side_for_token(token: str, recipe_name: str) -> str:
    role = token.split(".", 1)[0]

    if recipe_name == "AT":
        return "left" if role == "A" else "right"

    if recipe_name == "D_T":
        return "left" if role == "D" else "right"

    if recipe_name == "AG_IX":
        return "left" if role == "A" else "right"

    if recipe_name == "GU_XXVIII":
        return "left" if role == "G" else "right"

    # GC-like recipes use G-like geometry for recipe role 1 and
    # C-like geometry for recipe role 2, even for B:S, Z:P, K:X, etc.
    return "left" if role == "G" else "right"


def _atom_coord(residue, pair_residue: PairResidue, canonical: str, recipe_name: str):
    name = atom_name(pair_residue, canonical, recipe_name)
    return np.asarray(residue[name].coord, dtype=float)


def _angle(a, b, c) -> float:
    v1 = a - b
    v2 = c - b
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return 180.0
    x = np.clip(np.dot(v1, v2) / denom, -1.0, 1.0)
    return math.degrees(math.acos(x))


def _plane(residue, pair_residue: PairResidue, canonical_names, recipe_name):
    points = []

    for canonical in canonical_names:
        try:
            points.append(
                _atom_coord(residue, pair_residue, canonical, recipe_name)
            )
        except (KeyError, ValueError):
            continue

    if len(points) < 3:
        return None

    points = np.asarray(points)
    center = points.mean(axis=0)
    _, _, vh = np.linalg.svd(points - center)
    normal = vh[-1]
    normal /= np.linalg.norm(normal)

    return center, normal


def _plane_angle(plane1, plane2) -> float:
    x = abs(float(np.dot(plane1[1], plane2[1])))
    return math.degrees(math.acos(np.clip(x, -1.0, 1.0)))


def _candidate_score(
    r1,
    r2,
    p1: PairResidue,
    p2: PairResidue,
    recipe,
    *,
    allow_bad_geometry: bool = False,
) -> Candidate | None:
    roles = _role_residues(p1, p2, recipe)
    if roles is None:
        return None

    left_p, right_p = roles
    left_r, right_r = (r1, r2) if roles == (p1, p2) else (r2, r1)

    angles, bonds, left_plane_atoms, right_plane_atoms = _template(recipe.name)

    bond_distances = []

    for a1, a2, ideal, _sigma in bonds:
        role1, atom1 = a1.split(".", 1)
        role2, atom2 = a2.split(".", 1)

        side1 = _side_for_token(a1, recipe.name)
        side2 = _side_for_token(a2, recipe.name)

        try:
            c1 = _atom_coord(
                left_r if side1 == "left" else right_r,
                left_p if side1 == "left" else right_p,
                atom1,
                recipe.name,
            )
            c2 = _atom_coord(
                left_r if side2 == "left" else right_r,
                left_p if side2 == "left" else right_p,
                atom2,
                recipe.name,
            )
        except (KeyError, ValueError):
            return None

        bond_distances.append(float(np.linalg.norm(c1 - c2)))

    angle_errors = []

    for a1, a2, a3, ideal, _sigma in angles:
        coords = []

        for token in (a1, a2, a3):
            atom = token.split(".", 1)[1]
            side = _side_for_token(token, recipe.name)

            try:
                coords.append(
                    _atom_coord(
                        left_r if side == "left" else right_r,
                        left_p if side == "left" else right_p,
                        atom,
                        recipe.name,
                    )
                )
            except (KeyError, ValueError):
                return None

        angle_errors.append(abs(_angle(*coords) - ideal))

    left_plane = _plane(
        left_r, left_p, left_plane_atoms, recipe.name
    )
    right_plane = _plane(
        right_r, right_p, right_plane_atoms, recipe.name
    )

    if left_plane is None or right_plane is None:
        return None

    plane_angle = _plane_angle(left_plane, right_plane)
    center_distance = float(
        np.linalg.norm(left_plane[0] - right_plane[0])
    )

    if not allow_bad_geometry:
        if center_distance > SEARCH_RADIUS or plane_angle > MAX_PLANE_ANGLE:
            return None

    bond_score = sum(
        math.exp(-0.5 * ((d - ideal) / BOND_SOFT_SIGMA) ** 2)
        for d, (_, _, ideal, _) in zip(bond_distances, bonds)
    ) / len(bonds)

    angle_score = (
        sum(math.exp(-0.5 * (error / 12.0) ** 2) for error in angle_errors)
        / len(angle_errors)
        if angle_errors
        else 1.0
    )

    plane_score = math.exp(
        -0.5 * (plane_angle / 12.0) ** 2
    )

    center_score = math.exp(
        -0.5 * ((center_distance - 4.5) / 1.5) ** 2
    )

    score = 100.0 * (
        0.50 * bond_score
        + 0.25 * angle_score
        + 0.15 * plane_score
        + 0.10 * center_score
    )

    # Non-canonical candidates are deliberately de-emphasized.
    if recipe.noncanonical:
        score -= 10.0

    return Candidate(
        first=p1,
        second=p2,
        recipe=recipe.name,
        score=score,
        bond_distances=tuple(bond_distances),
        angle_errors=tuple(angle_errors),
        plane_angle=plane_angle,
        noncanonical=recipe.noncanonical,
    )


def _compatible_candidates(residues, *, allow_noncanonical: bool = False):
    for i, (r1, p1) in enumerate(residues):
        for r2, p2 in residues[i + 1:]:
            # Adjacent residues on the same chain are stacking neighbors,
            # not candidates for base pairing.
            if (
                r1.get_parent().id == r2.get_parent().id
                and abs(r1.id[1] - r2.id[1]) <= 1
            ):
                continue

            recipe = recipe_for(p1.base_class, p2.base_class)
            if recipe is None:
                continue
            if recipe.noncanonical and not allow_noncanonical:
                continue

            candidate = _candidate_score(
                r1, r2, p1, p2, recipe
            )

            if candidate is not None:
                yield candidate


def _run_key(candidate: Candidate):
    try:
        a = int(candidate.first.resid)
        b = int(candidate.second.resid)
    except ValueError:
        return None
    return (
        candidate.first.chain,
        candidate.second.chain,
        a,
        b,
    )


def _can_extend(left: Candidate, right: Candidate) -> bool:
    """Whether right continues left on antiparallel tracks."""
    if (
        left.first.chain != right.first.chain
        or left.second.chain != right.second.chain
    ):
        return False
    try:
        da = int(right.first.resid) - int(left.first.resid)
        db = int(right.second.resid) - int(left.second.resid)
    except ValueError:
        return False
    return da == 1 and db == -1


def _build_runs(candidates: list[Candidate]) -> list[list[Candidate]]:
    """Build contiguous antiparallel runs from plausible candidates."""
    by_key = {}
    for c in candidates:
        key = _run_key(c)
        if key is not None:
            by_key[key] = c

    runs = []
    used = set()

    for c in sorted(candidates, key=lambda x: x.score, reverse=True):
        key = _run_key(c)
        if key is None or key in used:
            continue

        run = [c]
        used.add(key)

        # Extend forward.
        current = c
        while True:
            try:
                a = int(current.first.resid)
                b = int(current.second.resid)
            except ValueError:
                break

            next_key = (
                current.first.chain,
                current.second.chain,
                a + 1,
                b - 1,
            )
            nxt = by_key.get(next_key)
            if nxt is None or next_key in used:
                break
            run.append(nxt)
            used.add(next_key)
            current = nxt

        # Extend backward.
        current = c
        backwards = []
        while True:
            try:
                a = int(current.first.resid)
                b = int(current.second.resid)
            except ValueError:
                break

            prev_key = (
                current.first.chain,
                current.second.chain,
                a - 1,
                b + 1,
            )
            prv = by_key.get(prev_key)
            if prv is None or prev_key in used:
                break
            backwards.append(prv)
            used.add(prev_key)
            current = prv

        run = list(reversed(backwards)) + run
        runs.append(run)

    return runs


def _run_score(run: list[Candidate]) -> float:
    if not run:
        return 0.0

    # Continuous runs are much stronger evidence than isolated pairs.
    score = sum(c.score for c in run) / len(run)
    score += 10.0 * max(0, len(run) - 1)
    score += 5.0 * max(0, len(run) - 2)
    return score


def _infer_run_extensions(
    run: list[Candidate],
    residue_lookup: dict[tuple[str, int], tuple[object, PairResidue]],
    *,
    allow_noncanonical: bool = False,
    metal_bound: frozenset = frozenset(),
):
    """
    Iteratively extend an established antiparallel run toward both ends.

    A newly accepted extension becomes the new edge, so a weak pair can be
    rescued by a strong run and the run can continue outward until a genuine
    hard stop or model boundary is reached.
    """
    additions = []
    messages = []

    if not run:
        return additions, messages

    # Work on copies so the caller's original run remains unchanged.
    working_run = list(run)

    def try_extension(edge: Candidate, da: int, db: int, direction: str):
        try:
            a = int(edge.first.resid) + da
            b = int(edge.second.resid) + db
        except ValueError:
            return None, True

        key1 = (edge.first.chain, a)
        key2 = (edge.second.chain, b)

        if (edge.first.chain, str(a)) in metal_bound or (edge.second.chain, str(b)) in metal_bound:
            messages.append(
                f"RUN EDGE: {edge.first.chain}{a} <-> "
                f"{edge.second.chain}{b}: metal-mediated base; run stops"
            )
            return None, True

        if key1 not in residue_lookup or key2 not in residue_lookup:
            messages.append(
                f"RUN EDGE: {edge.first.chain}{a} <-> "
                f"{edge.second.chain}{b} is outside the model"
            )
            return None, True

        r1, p1 = residue_lookup[key1]
        r2, p2 = residue_lookup[key2]

        recipe = recipe_for(p1.base_class, p2.base_class)
        if recipe is None or (recipe.noncanonical and not allow_noncanonical):
            messages.append(
                f"RUN EDGE: {edge.first.chain}{a} <-> "
                f"{edge.second.chain}{b}: no configured recipe"
            )
            return None, True

        candidate = _candidate_score(
            r1,
            r2,
            p1,
            p2,
            recipe,
            allow_bad_geometry=True,
        )

        if candidate is None:
            messages.append(
                f"RUN EDGE: {edge.first.chain}{a} <-> "
                f"{edge.second.chain}{b}: geometry could not be checked"
            )
            return None, True

        # Run context permits much worse individual geometry than a seed
        # candidate, but we still require something physically plausible.
        if candidate.score >= 10.0:
            accepted = Candidate(
                candidate.first,
                candidate.second,
                candidate.recipe,
                candidate.score + 20.0,
                candidate.bond_distances,
                candidate.angle_errors,
                candidate.plane_angle,
                candidate.noncanonical,
            )
            messages.append(
                f"RUN EDGE: {edge.first.chain}{a} <-> "
                f"{edge.second.chain}{b} score={candidate.score:.1f} "
                f"accepted as {direction} extension"
            )
            return accepted, False

        messages.append(
            f"RUN EDGE: {edge.first.chain}{a} <-> "
            f"{edge.second.chain}{b} score={candidate.score:.1f} "
            f"hard stop"
        )
        return None, True

    # Repeatedly extend the forward edge.
    forward = working_run[-1]
    while True:
        candidate, stop = try_extension(
            forward, +1, -1, "forward"
        )
        if stop:
            break
        additions.append(candidate)
        working_run.append(candidate)
        forward = candidate

    # Repeatedly extend the backward edge.
    backward = working_run[0]
    while True:
        candidate, stop = try_extension(
            backward, -1, +1, "backward"
        )
        if stop:
            break
        additions.append(candidate)
        working_run.insert(0, candidate)
        backward = candidate

    return additions, messages

def _infer_run_gaps(
    run: list[Candidate],
    residue_lookup: dict[tuple[str, int], tuple[object, PairResidue]],
    *,
    allow_noncanonical: bool = False,
    metal_bound: frozenset = frozenset(),
):
    """Check one-residue holes inside an otherwise continuous run."""
    additions = []
    messages = []

    if len(run) < 2:
        return additions, messages

    ordered = sorted(run, key=lambda c: int(c.first.resid))

    for left, right in zip(ordered, ordered[1:]):
        try:
            da = int(right.first.resid) - int(left.first.resid)
            db = int(right.second.resid) - int(left.second.resid)
        except ValueError:
            continue

        if da != 2 or db != -2:
            continue

        a = int(left.first.resid) + 1
        b = int(left.second.resid) - 1

        key1 = (left.first.chain, a)
        key2 = (left.second.chain, b)

        if (left.first.chain, str(a)) in metal_bound or (left.second.chain, str(b)) in metal_bound:
            messages.append(
                f"RUN GAP: {left.first.chain}{a} <-> "
                f"{left.second.chain}{b}: metal-mediated base; left unpaired"
            )
            continue

        if key1 not in residue_lookup or key2 not in residue_lookup:
            continue

        r1, p1 = residue_lookup[key1]
        r2, p2 = residue_lookup[key2]

        recipe = recipe_for(p1.base_class, p2.base_class)
        if recipe is None or (recipe.noncanonical and not allow_noncanonical):
            messages.append(
                f"RUN GAP: {left.first.chain}{a} <-> "
                f"{left.second.chain}{b}: no configured recipe"
            )
            continue

        candidate = _candidate_score(r1, r2, p1, p2, recipe)

        if candidate is None:
            messages.append(
                f"RUN GAP: {left.first.chain}{a} <-> "
                f"{left.second.chain}{b}: geometry could not be checked"
            )
            continue

        if candidate.score >= 30.0:
            additions.append(
                Candidate(
                    candidate.first,
                    candidate.second,
                    candidate.recipe,
                    candidate.score + 25.0,
                    candidate.bond_distances,
                    candidate.angle_errors,
                    candidate.plane_angle,
                )
            )
            messages.append(
                f"RUN GAP: {left.first.chain}{a} <-> "
                f"{left.second.chain}{b} score={candidate.score:.1f} "
                f"accepted"
            )
        else:
            messages.append(
                f"RUN GAP: {left.first.chain}{a} <-> "
                f"{left.second.chain}{b} score={candidate.score:.1f} "
                f"hard stop"
            )

    return additions, messages


def _select_runs(
    candidates: list[Candidate],
    residue_lookup: dict[tuple[str, int], tuple[object, PairResidue]],
    *,
    allow_noncanonical: bool = False,
    metal_bound: frozenset = frozenset(),
):
    """Select strong antiparallel runs, then fill gaps and extend edges."""
    # Only reasonably plausible candidates seed runs. We keep the threshold
    # lower than the final isolated-pair threshold because continuity matters.
    seeds = [
        c for c in candidates
        if c.score >= (55.0 if c.noncanonical else 45.0)
    ]
    runs = _build_runs(seeds)

    runs.sort(key=lambda r: (_run_score(r), len(r)), reverse=True)

    selected = []
    used = set()
    messages = []

    for run in runs:
        run_additions, run_messages = _infer_run_gaps(
            run, residue_lookup, allow_noncanonical=allow_noncanonical,
            metal_bound=metal_bound,
        )
        edge_additions, edge_messages = _infer_run_extensions(
            run, residue_lookup, allow_noncanonical=allow_noncanonical,
            metal_bound=metal_bound,
        )

        messages.extend(run_messages)
        messages.extend(edge_messages)

        proposed = list(run)
        proposed.extend(run_additions)
        proposed.extend(edge_additions)

        # Prefer the best available candidate at each residue position.
        proposed.sort(key=lambda c: c.score, reverse=True)

        accepted = []
        local_used = set()

        for candidate in proposed:
            a = (candidate.first.chain, candidate.first.resid)
            b = (candidate.second.chain, candidate.second.resid)

            if a in local_used or b in local_used:
                continue
            if a in used or b in used:
                continue

            accepted.append(candidate)
            local_used.update((a, b))

        # A run needs at least two internally compatible pairs to establish
        # enough structural context for its rescued additions.
        if len(accepted) < 2:
            continue

        selected.extend(accepted)
        used.update(local_used)

    # Finally allow isolated, high-confidence pairs not belonging to a run.
    for candidate in sorted(candidates, key=lambda c: c.score, reverse=True):
        a = (candidate.first.chain, candidate.first.resid)
        b = (candidate.second.chain, candidate.second.resid)

        if candidate.score < (55.0 if candidate.noncanonical else MIN_SCORE):
            continue
        if a in used or b in used:
            continue

        selected.append(candidate)
        used.update((a, b))

    selected.sort(
        key=lambda c: (
            c.first.chain,
            int(c.first.resid),
            c.second.chain,
            int(c.second.resid),
        )
    )

    return selected, messages


def guess_pairs(
    pdb_filename: str | Path,
    *,
    allow_noncanonical: bool = False,
    metal_pairs: bool = True,
    metal_distance: float = DEFAULT_MAX_DISTANCE,
    metal_elements=None,
) -> tuple[list[Candidate], list[str]]:
    """Guess base pairs; see :func:`guess_pairs_detailed`."""
    selected, messages, _ = guess_pairs_detailed(
        pdb_filename,
        allow_noncanonical=allow_noncanonical,
        metal_pairs=metal_pairs,
        metal_distance=metal_distance,
        metal_elements=metal_elements,
    )
    return selected, messages


def guess_pairs_detailed(
    pdb_filename: str | Path,
    *,
    allow_noncanonical: bool = False,
    metal_pairs: bool = True,
    metal_distance: float = DEFAULT_MAX_DISTANCE,
    metal_elements=None,
) -> tuple[list[Candidate], list[str], list[MetalMediatedPair]]:
    """Guess base pairs and report metal-mediated pairs separately.

    Bases in metal-mediated pairs (``restraints.metal_pairs``) are never
    offered as candidates, never used to extend a run and never fill a gap,
    so they cannot be guessed as a canonical pair with anything.
    """
    found: list[MetalMediatedPair] = []
    metal_messages: list[str] = []
    if metal_pairs:
        found, metal_messages = find_metal_mediated_pairs(
            pdb_filename, max_distance=metal_distance, metals=metal_elements
        )
    metal_bound = frozenset(metal_bound_residues(found))

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(
        "NARestraints", str(pdb_filename)
    )
    records = load_residue_records()

    residues = []

    for model in first_model(structure):
        for chain in model:
            for residue in chain.get_residues():
                if (chain.id, str(residue.id[1])) in metal_bound:
                    continue
                try:
                    pair_residue = pair_residue_from_pdb(
                        records, residue
                    )
                except ValueError:
                    continue

                residues.append((residue, pair_residue))

    candidates = list(
        _compatible_candidates(residues, allow_noncanonical=allow_noncanonical)
    )

    residue_lookup = {
        (p.chain, int(p.resid)): (r, p)
        for r, p in residues
        if p.resid.isdigit()
    }

    selected, messages = _select_runs(
        candidates, residue_lookup, allow_noncanonical=allow_noncanonical,
        metal_bound=metal_bound,
    )
    return selected, metal_messages + messages, found


def write_guess(
    filename: str | Path,
    candidates: Iterable[Candidate],
    metal_pairs: Iterable[MetalMediatedPair] = (),
) -> None:
    """Write guessed pairs; metal-mediated pairs are written as comments.

    The commented blocks are ignored by ``narestraints``, so these pairs are
    recorded for the reader without being restrained.
    """
    lines = []

    for candidate in candidates:
        lines.append(
            f"{candidate.first.chain} {candidate.first.resid}"
        )
        lines.append(
            f"{candidate.second.chain} {candidate.second.resid}"
        )
        if candidate.noncanonical:
            lines.append("# non-canonical: true")
            if any(d > 3.5 for d in candidate.bond_distances):
                lines.append(
                    "# WARNING: non-canonical H-bond geometry is currently poor; "
                    "this may be expected before refinement"
            )
        lines.append("")

    for pair in metal_pairs:
        first, second = pair.residues
        lines.append(
            "# Metal-mediated pair: not a canonical pair and not restrained by NARestraints"
        )
        lines.append(f"# {first[0].strip()} {first[1]}")
        lines.append(f"# {second[0].strip()} {second[1]}")
        lines.append(f"# {pair.describe()}")
        lines.append("")

    Path(filename).write_text(
        "\n".join(lines).rstrip() + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Guess nucleic-acid base pairs from low-resolution "
            "PDB geometry. This is independent of restraint generation."
        )
    )
    parser.add_argument("pdb", help="Input PDB file")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output .txt filename (default: <PDB stem>-guessed.txt)",
    )
    parser.add_argument(
        "--expected-pairs",
        type=int,
        default=None,
        help="Expected number of base pairs; mismatch produces a warning",
    )
    parser.add_argument(
        "--non-canonical",
        action="store_true",
        help="Explicitly allow non-canonical pair recipes (AG IX and GU XXVIII)",
    )
    from .cli import (
        add_metal_pair_arguments,
        add_readyset_arguments,
        metal_elements_from_args,
        prepared_model_from_args,
    )

    add_metal_pair_arguments(parser)
    add_readyset_arguments(parser)

    args = parser.parse_args()
    metal_elements = metal_elements_from_args(parser, args)

    pdb_path = Path(args.pdb)
    output = (
        args.output
        or pdb_path.with_name(
            pdb_path.stem + "-guessed.txt"
        )
    )

    model = prepared_model_from_args(parser, args)
    fixed_pdb = validate_and_fix_pdb(model)
    candidates, messages, metal_pairs = guess_pairs_detailed(
        fixed_pdb,
        allow_noncanonical=args.non_canonical,
        metal_pairs=not args.no_metal_pairs,
        metal_distance=args.metal_distance,
        metal_elements=metal_elements,
    )

    print(f"BASE PAIRS FOUND: {len(candidates)}")
    for pair in metal_pairs:
        print(f"METAL-MEDIATED PAIR (not guessed as canonical): {pair.describe()}")

    if args.expected_pairs is not None:
        if len(candidates) != args.expected_pairs:
            print(
                f"WARNING: expected {args.expected_pairs} "
                f"base pairs, found {len(candidates)}"
            )
        else:
            print(
                f"EXPECTED COUNT MATCHES: {args.expected_pairs}"
            )

    for message in messages:
        print(message)

    for candidate in candidates:
        print(
            f"  {candidate.first.chain}{candidate.first.resid} "
            f"<-> {candidate.second.chain}{candidate.second.resid} "
            f"{candidate.recipe:10s} "
            f"score={candidate.score:5.1f} "
            f"bonds="
            + ",".join(
                f"{distance:.2f}"
                for distance in candidate.bond_distances
            )
            + f" plane={candidate.plane_angle:.1f}°"
            + (" non-canonical=true" if candidate.noncanonical else "")
            + (
                " WARNING: poor current H-bond geometry; check after refinement"
                if candidate.noncanonical
                and any(d > 3.5 for d in candidate.bond_distances)
                else ""
            )
        )

    write_guess(output, candidates, metal_pairs)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
