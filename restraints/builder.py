from __future__ import annotations

from pathlib import Path
from typing import Iterable

from Bio.PDB import PDBParser

from .anomalous import (
    DEFAULT_MIN_ATOMIC_NUMBER,
    HeavyAtomSite,
    anomalous_phil_block,
    describe_sites,
    find_heavy_atoms,
)
from .base_pairs import BasePairStretch
from .metal_pairs import (
    DEFAULT_MAX_DISTANCE,
    find_metal_mediated_pairs,
    metal_bound_residues,
)
from .phenix import PairResidue, generate_pair_restraints, write_phil
from .residue_library import find_residue, load_residue_records
from .stacking import StackingResidue, generate_chain_stacking_restraints


_PAIR_CATEGORIES = {"A", "T", "G", "C", "D", "B", "S", "Z", "P", "K", "X", "I"}


def _record_score(record: dict, residue) -> int:
    """Prefer the residue-library mapping that actually matches the PDB atoms."""
    score = 0
    for key, value in record.items():
        if key in {
            "Ligand code", "Name", "Abbreviation", "Base Analog", "Phosphate",
            "Sugar Type", "SugarType", "Saenger", "Notes", "Ref", "Entries",
            "Function", "Source sheet", "Glycosidic atom",
        }:
            continue
        if value is not None and value != "" and residue.has_id(str(value).strip()):
            score += 1
    return score


def record_for_residue(records: list[dict], residue) -> dict:
    code = residue.get_resname().strip()
    matches = find_residue(records, code)
    if not matches:
        raise ValueError(
            f"Residue {code!r} is not present in the residue library."
        )
    return max(matches, key=lambda record: _record_score(record, residue))


def pair_residue_from_pdb(records: list[dict], residue) -> PairResidue:
    record = record_for_residue(records, residue)
    base_class = record.get("Base Analog")

    if base_class not in _PAIR_CATEGORIES:
        raise ValueError(
            f"{residue.get_resname().strip()} maps to unsupported Base Analog "
            f"{base_class!r}. Add it to the recipe/category configuration first."
        )

    atoms = {}
    for key, value in record.items():
        if key in {
            "Ligand code", "Name", "Abbreviation", "Base Analog", "Phosphate",
            "Sugar Type", "SugarType", "Saenger", "Notes", "Ref", "Entries",
            "Function", "Source sheet", "Glycosidic atom",
        }:
            continue
        if value is not None and value != "" and residue.has_id(str(value).strip()):
            atoms[key] = str(value).strip()

    return PairResidue(
        chain=residue.get_parent().get_id(),
        resid=str(residue.id[1]),
        atoms=atoms,
        base_class=base_class,
    )


def _get_residue(structure, chain: str, resid: str):
    model = next(structure.get_models())
    chain_obj = model[chain]
    # PDB insertion codes are retained if present; current input format is integer resid.
    for residue in chain_obj.get_residues():
        if str(residue.id[1]) == str(resid):
            return residue
    raise ValueError(f"Could not find {chain}:{resid} in the PDB.")


def _stacking_residue(records: list[dict], residue) -> StackingResidue:
    """Build a stacking record without making unknown chemistry fatal.

    Known residues receive the full workbook mapping needed for explicit
    modified-base stacking. If a residue is not yet mappable, retain its
    identity so the stacking layer can fall back to a generic Phenix
    ``stacking_pair`` and emit a warning rather than dropping the stack.
    """
    mapped = None
    try:
        mapped = pair_residue_from_pdb(records, residue)
    except ValueError:
        pass

    return StackingResidue(
        chain=residue.get_parent().get_id(),
        resid=str(residue.id[1]),
        residue_name=residue.get_resname().strip(),
        mapped=mapped,
    )


def _is_nucleotide_like(residue, stacking: StackingResidue) -> bool:
    """Whether a residue can take part in base stacking.

    Workbook-mapped residues and any residue with a C1' (or C1*) atom count.
    Ions, waters and other ligands do not, so they are never stacked onto a
    neighbouring base just because they follow it in the chain.
    """
    return (
        stacking.mapped is not None
        or residue.has_id("C1'")
        or residue.has_id("C1*")
    )


def _residue_key(chain: str, resid) -> tuple[str, str]:
    return (chain, str(int(str(resid))))


def _split_at(residues, excluded: set) -> list[list[StackingResidue]]:
    """Split a chain into stretches that skip the excluded residues.

    Splitting (rather than filtering) ensures no new stack is created across
    an excluded residue between its two neighbours.
    """
    stretches: list[list[StackingResidue]] = [[]]
    for residue in residues:
        if _residue_key(residue.chain, residue.resid) in excluded:
            stretches.append([])
        else:
            stretches[-1].append(residue)
    return [stretch for stretch in stretches if len(stretch) >= 2]


def _summarise_names(names: list[str]) -> str:
    counts: dict[str, int] = {}
    for name in names:
        counts[name] = counts.get(name, 0) + 1
    return ", ".join(f"{name} x{count}" for name, count in counts.items())


def build_phil_from_pdb(
    pdb_filename: str | Path,
    stretches: Iterable[BasePairStretch],
    output_filename: str | Path,
    *,
    parallels: bool = True,
    planes: bool = True,
    include_stacking: bool = True,
    anomalous: bool = False,
    anomalous_elements: Iterable[str] | None = None,
    anomalous_strategy: bool = True,
    anomalous_min_atomic_number: int = DEFAULT_MIN_ATOMIC_NUMBER,
    metal_pairs: bool = True,
    metal_distance: float = DEFAULT_MAX_DISTANCE,
    metal_elements: Iterable[str] | None = None,
    keep_metal_pair_stacking: bool = False,
) -> list[HeavyAtomSite]:
    """Write the restraint PHIL and return any anomalous scatterers written.

    With ``anomalous=True`` (or an explicit ``anomalous_elements`` list),
    heavy atoms in HETATM records are written as phenix.refine anomalous
    scatterer groups; see ``restraints.anomalous``. ``anomalous_strategy``
    also adds ``group_anomalous`` to the refinement strategy, without which
    phenix.refine does not refine the groups.

    With ``metal_pairs=True`` (the default), metal-mediated base pairs are
    found from the coordinates (see ``restraints.metal_pairs``). Any pair-file
    pair involving one of their bases is skipped, and so is stacking of those
    bases unless ``keep_metal_pair_stacking`` is set.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("NARestraints", str(pdb_filename))
    records = load_residue_records()

    bound: dict = {}
    if metal_pairs:
        found_pairs, metal_messages = find_metal_mediated_pairs(
            pdb_filename, max_distance=metal_distance, metals=metal_elements
        )
        for message in metal_messages:
            print(f"NOTE: {message}")
        for found in found_pairs:
            print(f"Metal-mediated pair (not restrained): {found.describe()}")
        bound = metal_bound_residues(found_pairs)

    pair_blocks: list[str] = []

    for stretch in stretches:
        for pair in stretch.pairs():
            key1 = _residue_key(pair.base1.chain, pair.base1.resid)
            key2 = _residue_key(pair.base2.chain, pair.base2.resid)
            metal_pair = bound.get(key1) or bound.get(key2)
            if metal_pair is not None:
                listed = f"{pair.base1.chain}:{pair.base1.resid} <-> {pair.base2.chain}:{pair.base2.resid}"
                if {key1, key2} == set(metal_pair.residues):
                    print(
                        f"Skipping base-pair restraints for {listed}: "
                        f"metal-mediated via {metal_pair.metal_label}."
                    )
                else:
                    print(
                        f"WARNING: Skipping base-pair restraints for {listed}: "
                        f"a base of this pair is in the metal-mediated pair "
                        f"{metal_pair.describe()}. Check the pair file."
                    )
                continue

            r1 = _get_residue(structure, pair.base1.chain, pair.base1.resid)
            r2 = _get_residue(structure, pair.base2.chain, pair.base2.resid)

            p1 = pair_residue_from_pdb(records, r1)
            p2 = pair_residue_from_pdb(records, r2)

            block = generate_pair_restraints(
                p1,
                p2,
                parallels=parallels,
                planes=planes,
            )
            if block is None:
                print(
                    "WARNING: No restraint recipe for pair "
                    f"{p1.chain}:{p1.resid} ({p1.base_class}) <-> "
                    f"{p2.chain}:{p2.resid} ({p2.base_class}); "
                    "skipping pair-specific restraints."
                )
                continue
            pair_blocks.append(block)

    stacking_block = ""
    if include_stacking:
        # Match the MATLAB behavior: stacking is generated for consecutive
        # residues in each PDB chain, independently of selected base pairs.
        # Canonical stacks remain Phenix stacking_pair records. Any stack
        # touching a mapped modified residue becomes an explicit parallelity
        # geometry edit so Phenix does not need to infer a plane from the
        # modified residue name.
        #
        # Only nucleotide-like residues are stacked: ions, waters and other
        # ligands that share a chain with the nucleic acid are skipped. Bases
        # in metal-mediated pairs are not stacked unless explicitly kept.
        excluded = set() if keep_metal_pair_stacking else set(bound)
        for model in structure:
            for chain in model:
                stacking_residues = []
                skipped_names = []
                for residue in chain.get_residues():
                    stacking = _stacking_residue(records, residue)
                    if _is_nucleotide_like(residue, stacking):
                        stacking_residues.append(stacking)
                    else:
                        skipped_names.append(stacking.residue_name)
                if skipped_names:
                    print(
                        f"Chain {chain.id}: no stacking for non-nucleotide residues "
                        f"({_summarise_names(skipped_names)})."
                    )

                skipped_bound = [
                    f"{residue.chain}:{residue.resid}"
                    for residue in stacking_residues
                    if _residue_key(residue.chain, residue.resid) in excluded
                ]
                if skipped_bound:
                    print(
                        f"Chain {chain.id}: no stacking for metal-mediated base(s) "
                        f"{', '.join(skipped_bound)}."
                    )

                for stretch in _split_at(stacking_residues, excluded):
                    generic, manual, warnings = generate_chain_stacking_restraints(
                        stretch
                    )
                    if generic:
                        stacking_block += (
                            "\n\n" if stacking_block else ""
                        ) + generic
                    pair_blocks.extend(manual)
                    for warning in warnings:
                        print(f"WARNING: {warning}")

    anomalous_sites: list[HeavyAtomSite] = []
    anomalous_block = ""
    if anomalous or anomalous_elements is not None:
        anomalous_sites, anomalous_warnings = find_heavy_atoms(
            pdb_filename,
            elements=anomalous_elements,
            min_atomic_number=anomalous_min_atomic_number,
        )
        for warning in anomalous_warnings:
            print(f"WARNING: {warning}")
        if anomalous_sites:
            for line in describe_sites(anomalous_sites):
                print(line)
            anomalous_block = anomalous_phil_block(
                anomalous_sites, include_strategy=anomalous_strategy
            )
        else:
            print(
                "WARNING: --anomalous found no heavy atoms in HETATM records; "
                "no anomalous scatterer groups were written."
            )

    write_phil(output_filename, pair_blocks, stacking_block, anomalous_block)
    return anomalous_sites
