
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from .recipe_library import PairRecipe, recipe_for


AT_ANGLE_RESTRAINTS = (
    # atom 1, atom 2, atom 3, ideal, sigma
    ("A.C2", "A.N1", "T.N3", 116.2, 3.46),
    ("A.C6", "A.N6", "T.O4", 115.6, 8.34),
)

GC_ANGLE_RESTRAINTS = (
    ("G.C2", "G.N1", "C.N3", 116.2, 3.46),
    ("G.C2", "G.N2", "C.O2", 122.2, 2.88),
    ("C.C4", "C.N4", "G.O6", 117.3, 2.86),
    ("C.C2", "C.O2", "G.N2", 120.7, 2.20),
    ("C.C2", "C.N3", "G.N1", 115.8, 2.88),
    ("G.C6", "G.O6", "C.N4", 122.8, 3.00),
)

AT_BONDS = (
    ("T.N3", "A.N1", 2.8, 0.2),
    ("T.O4", "A.N6", 2.8, 0.2),
)

GC_BONDS = (
    ("G.N2", "C.O2", 2.8, 0.2),
    ("G.N1", "C.N3", 2.8, 0.2),
    ("G.O6", "C.N4", 2.8, 0.2),
)

DT_ANGLE_RESTRAINTS = (
    ("D.C2", "D.N1", "T.N3", 116.2, 3.46),
    ("D.C6", "D.N6", "T.O4", 115.6, 8.34),
    ("D.C2", "D.N2", "T.O2", 122.2, 2.88),
    ("T.C2", "T.O2", "D.N2", 120.7, 2.20),
)

DT_BONDS = (
    ("D.N1", "T.N3", 2.8, 0.2),
    ("D.N6", "T.O4", 2.8, 0.2),
    ("D.N2", "T.O2", 2.8, 0.2),
)

# Non-canonical pair contacts. These deliberately use the experimentally
# established hydrogen-bond distances only; recipe-specific angle targets
# should be added later from validated structural statistics.
AG_IX_ANGLES = ()
AG_IX_BONDS = (
    ("A.N6", "G.O6", 2.8, 0.2),
    ("A.N1", "G.N1", 2.8, 0.2),
)

GU_XXVIII_ANGLES = ()
GU_XXVIII_BONDS = (
    ("G.O6", "T.N3", 2.8, 0.2),
    ("G.N1", "T.O2", 2.8, 0.2),
)

AT_PARALLEL_T = ("C2", "O2", "N1", "N3", "C4", "O4", "C5", "C7", "C6")
AT_PARALLEL_A = ("C2", "N1", "C6", "N6", "C5", "C4", "N3", "N9", "C8", "N7")

GC_PARALLEL_G = ("C2", "N2", "C6", "O6", "N1", "C5", "C4", "N3", "N9", "C8", "N7")
GC_PARALLEL_C = ("C2", "O2", "N1", "N3", "C4", "N4", "C5", "C6")


@dataclass(frozen=True)
class PairResidue:
    chain: str
    resid: str
    atoms: dict[str, str]
    base_class: str


# Recipe-role aliases resolve a canonical atom role requested by a recipe to
# the actual atom role exposed by a residue category.  These are deliberately
# explicit: they describe how a residue category can fill a role in a recipe,
# rather than changing the underlying Ligands.xlsx mapping.
RECIPE_ROLE_ALIASES = {
    "GC": {
        "T": {"N4": "O4"},
    },
}


def atom_name(residue: PairResidue, canonical: str, recipe_name: str | None = None) -> str:
    candidates = [canonical]
    if recipe_name is not None:
        alias = RECIPE_ROLE_ALIASES.get(recipe_name, {}).get(residue.base_class, {})
        if canonical in alias:
            candidates.append(alias[canonical])

    for candidate in candidates:
        if candidate in residue.atoms:
            return residue.atoms[candidate]

    raise ValueError(
        f"{residue.chain}:{residue.resid} ({residue.base_class}) "
        f"has no mapped atom {canonical!r} for {recipe_name or 'the requested'} recipe"
    )


def atom_selection(
    residue: PairResidue, canonical: str, recipe_name: str | None = None
) -> str:
    return (
        f"chain {residue.chain} and resid {residue.resid} "
        f"and name {atom_name(residue, canonical, recipe_name)}"
    )


def _selection_for_atoms(
    residue: PairResidue, canonical_names: Sequence[str], recipe_name: str | None = None
) -> str:
    names = []
    for canonical in canonical_names:
        try:
            names.append(atom_name(residue, canonical, recipe_name))
        except ValueError:
            continue
    if not names:
        raise ValueError(
            f"No requested ring atoms are present for "
            f"{residue.chain}:{residue.resid}"
        )
    joined = " or ".join(f"name {name}" for name in names)
    return f"chain {residue.chain} and resid {residue.resid} and ({joined})"


def _angle_block(a1: str, a2: str, a3: str, ideal: float, sigma: float) -> str:
    return "\n".join(
        [
            "    angle {",
            "      action = add",
            f"      atom_selection_1 = {a1}",
            f"      atom_selection_2 = {a2}",
            f"      atom_selection_3 = {a3}",
            f"      angle_ideal = {ideal:g}",
            f"      sigma = {sigma:g}",
            "    }",
        ]
    )


def _bond_block(a1: str, a2: str, ideal: float, sigma: float) -> str:
    return "\n".join(
        [
            "    bond {",
            "      action=add",
            f"      atom_selection_1 = {a1}",
            f"      atom_selection_2 = {a2}",
            f"      distance_ideal = {ideal:g}",
            f"      sigma = {sigma:g}",
            "    }",
        ]
    )


def _parallelity_block(sel1: str, sel2: str) -> str:
    return "\n".join(
        [
            "    parallelity {",
            "      action = add",
            f"      atom_selection_1 = {sel1}",
            f"      atom_selection_2 = {sel2}",
            "      sigma = 0.027",
            "      target_angle_deg = 0",
            "    }",
        ]
    )


def _planarity_block(sel1: str, sel2: str) -> str:
    # Keep the MATLAB syntax: a single planarity group containing both bases.
    return "\n".join(
        [
            "    planarity {",
            "      action = add",
            f"      atom_selection = ({sel1}) or ({sel2})",
            "      sigma = 0.176",
            "    }",
        ]
    )


def _resolve_pair_orientation(
    first: PairResidue, second: PairResidue, recipe: PairRecipe
) -> tuple[PairResidue, PairResidue]:
    """Return residues in the role order required by the recipe."""
    if first.base_class in recipe.role1_categories and second.base_class in recipe.role2_categories:
        return first, second
    if second.base_class in recipe.role1_categories and first.base_class in recipe.role2_categories:
        return second, first
    raise ValueError(
        f"Pair categories {first.base_class}-{second.base_class} do not fit "
        f"the configured {recipe.name} recipe."
    )


def _bond_parameters(
    residue1: PairResidue,
    canonical1: str,
    residue2: PairResidue,
    canonical2: str,
    ideal: float,
    sigma: float,
    recipe_name: str,
    sulfur_pair: bool,
) -> tuple[float, float]:
    """Set distance/sigma for a bond in a sulfur-containing base pair.

    The sulfur contact itself is the strong anchor: 3.04 Å, sigma 0.2.
    The other contacts in the same pair keep their ideal distance but use
    sigma 0.4 Å so they can redistribute around the sulfur contact.
    """
    actual1 = atom_name(residue1, canonical1, recipe_name)
    actual2 = atom_name(residue2, canonical2, recipe_name)
    sulfur_contact = actual1.startswith("S") or actual2.startswith("S")

    if sulfur_pair:
        if sulfur_contact:
            return 3.04, 0.2
        return ideal, 0.4

    return ideal, sigma


def generate_pair_restraints(
    first: PairResidue,
    second: PairResidue,
    *,
    parallels: bool = True,
    planes: bool = True,
) -> str | None:
    """Generate restraints for one pair.

    Unsupported pairs return None. They are intentionally not fatal: the
    caller reports the missing recipe and continues generating the .phil file.
    """
    recipe = recipe_for(first.base_class, second.base_class)
    if recipe is None:
        return None

    left, right = _resolve_pair_orientation(first, second, recipe)
    lines: list[str] = []

    if recipe.name == "AT":
        refs = {"A": left, "T": right}
        for a1, a2, a3, ideal, sigma in AT_ANGLE_RESTRAINTS:
            lines.append(
                _angle_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    atom_selection(refs[a3[0]], a3[2:], recipe.name),
                    ideal,
                    sigma,
                )
            )
        for a1, a2, ideal, sigma in AT_BONDS:
            lines.append(
                _bond_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    *_bond_parameters(
                        refs[a1[0]], a1[2:],
                        refs[a2[0]], a2[2:],
                        ideal, sigma, recipe.name,
                        any(
                            atom_name(ref, canonical, recipe.name).startswith("S")
                            for ref in refs.values()
                            for canonical in ref.atoms
                        ),
                    ),
                )
            )
        if parallels:
            lines.append(
                _parallelity_block(
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                    _selection_for_atoms(left, AT_PARALLEL_A, recipe.name),
                )
            )
        if planes:
            lines.append(
                _planarity_block(
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                    _selection_for_atoms(left, AT_PARALLEL_A, recipe.name),
                )
            )
    elif recipe.name == "D_T":
        refs = {"D": left, "T": right}
        for a1, a2, a3, ideal, sigma in DT_ANGLE_RESTRAINTS:
            lines.append(
                _angle_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    atom_selection(refs[a3[0]], a3[2:], recipe.name),
                    ideal,
                    sigma,
                )
            )
        for a1, a2, ideal, sigma in DT_BONDS:
            lines.append(
                _bond_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    *_bond_parameters(
                        refs[a1[0]], a1[2:],
                        refs[a2[0]], a2[2:],
                        ideal, sigma, recipe.name,
                        any(
                            atom_name(ref, canonical, recipe.name).startswith("S")
                            for ref in refs.values()
                            for canonical in ref.atoms
                        ),
                    ),
                )
            )
        if parallels:
            lines.append(
                _parallelity_block(
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                )
            )
        if planes:
            lines.append(
                _planarity_block(
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                )
            )
    elif recipe.name == "AG_IX":
        refs = {"A": left, "G": right}
        for a1, a2, a3, ideal, sigma in AG_IX_ANGLES:
            lines.append(
                _angle_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    atom_selection(refs[a3[0]], a3[2:], recipe.name),
                    ideal, sigma,
                )
            )
        for a1, a2, ideal, sigma in AG_IX_BONDS:
            lines.append(
                _bond_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    *_bond_parameters(
                        refs[a1[0]], a1[2:],
                        refs[a2[0]], a2[2:],
                        ideal, sigma, recipe.name, False,
                    ),
                )
            )
        if parallels:
            lines.append(
                _parallelity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, GC_PARALLEL_G, recipe.name),
                )
            )
        if planes:
            lines.append(
                _planarity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, GC_PARALLEL_G, recipe.name),
                )
            )

    elif recipe.name == "GU_XXVIII":
        refs = {"G": left, "T": right}
        for a1, a2, a3, ideal, sigma in GU_XXVIII_ANGLES:
            lines.append(
                _angle_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    atom_selection(refs[a3[0]], a3[2:], recipe.name),
                    ideal, sigma,
                )
            )
        for a1, a2, ideal, sigma in GU_XXVIII_BONDS:
            lines.append(
                _bond_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    *_bond_parameters(
                        refs[a1[0]], a1[2:],
                        refs[a2[0]], a2[2:],
                        ideal, sigma, recipe.name, False,
                    ),
                )
            )
        if parallels:
            lines.append(
                _parallelity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                )
            )
        if planes:
            lines.append(
                _planarity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, AT_PARALLEL_T, recipe.name),
                )
            )

    else:
        refs = {"G": left, "C": right}
        for a1, a2, a3, ideal, sigma in GC_ANGLE_RESTRAINTS:
            lines.append(
                _angle_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    atom_selection(refs[a3[0]], a3[2:], recipe.name),
                    ideal,
                    sigma,
                )
            )
        for a1, a2, ideal, sigma in GC_BONDS:
            lines.append(
                _bond_block(
                    atom_selection(refs[a1[0]], a1[2:], recipe.name),
                    atom_selection(refs[a2[0]], a2[2:], recipe.name),
                    *_bond_parameters(
                        refs[a1[0]], a1[2:],
                        refs[a2[0]], a2[2:],
                        ideal, sigma, recipe.name,
                        any(
                            atom_name(ref, canonical, recipe.name).startswith("S")
                            for ref in refs.values()
                            for canonical in ref.atoms
                        ),
                    ),
                )
            )
        if parallels:
            lines.append(
                _parallelity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, GC_PARALLEL_C, recipe.name),
                )
            )
        if planes:
            lines.append(
                _planarity_block(
                    _selection_for_atoms(left, GC_PARALLEL_G, recipe.name),
                    _selection_for_atoms(right, GC_PARALLEL_C, recipe.name),
                )
            )

    return "\n\n".join(lines)

def generate_stacking_restraints(
    residue_sequence: Sequence[tuple[str, str]],) -> str:
    """Generate stacking restraints for consecutive residues in a chain."""
    blocks: list[str] = []

    for (chain1, resid1), (chain2, resid2) in zip(
        residue_sequence, residue_sequence[1:]
    ):
        if chain1 != chain2:
            continue

        blocks.append(
            "\n".join(
                [
                    "        stacking_pair {",
                    f"          base1 = chain {chain1} and resid {resid1}",
                    f"          base2 = chain {chain2} and resid {resid2}",
                    "        }",
                ]
            )
        )

    return "\n\n".join(blocks)

def write_phil(
    filename: str | Path,
    pair_blocks: Sequence[str],
    stacking_block: str = "",
    anomalous_block: str = "",
) -> None:
    """Write Phenix PHIL with custom geometry edits at the correct scope.

    IMPORTANT:
      refinement.geometry_restraints.edits is a sibling of
      refinement.pdb_interpretation, not a child of pdb_interpretation.

    ``anomalous_block`` is an optional ``refine { ... }`` scope (see
    ``restraints.anomalous.anomalous_phil_block``). It is written inside
    ``refinement`` after the geometry edits; when empty, the output is
    unchanged from earlier versions.
    """
    path = Path(filename)

    chunks = [
        "refinement {",
        "  pdb_interpretation {",
        "    secondary_structure {",
        "      nucleic_acid {",
    ]

    if stacking_block:
        chunks.append(stacking_block)

    chunks.extend([
        "      }",
        "      enabled = True",
        "    }",
        "  }",
        "",
        "  geometry_restraints.edits {",
    ])

    if pair_blocks:
        chunks.append("\n\n".join(pair_blocks))

    chunks.append("  }")

    if anomalous_block:
        chunks.extend(["", anomalous_block])

    chunks.append("}")

    path.write_text("\n".join(chunks) + "\n")
