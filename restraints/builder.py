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
) -> list[HeavyAtomSite]:
    """Write the restraint PHIL and return any anomalous scatterers written.

    With ``anomalous=True`` (or an explicit ``anomalous_elements`` list),
    heavy atoms in HETATM records are written as phenix.refine anomalous
    scatterer groups; see ``restraints.anomalous``. ``anomalous_strategy``
    also adds ``group_anomalous`` to the refinement strategy, without which
    phenix.refine does not refine the groups.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("NARestraints", str(pdb_filename))
    records = load_residue_records()

    pair_blocks: list[str] = []

    for stretch in stretches:
        for pair in stretch.pairs():
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
        for model in structure:
            for chain in model:
                stacking_residues = [
                    _stacking_residue(records, residue)
                    for residue in chain.get_residues()
                ]

                if len(stacking_residues) < 2:
                    continue

                generic, manual, warnings = generate_chain_stacking_restraints(
                    stacking_residues
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
