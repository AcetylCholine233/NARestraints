"""Find metal-mediated base pairs from coordinates.

A metal-mediated base pair (for example C-Ag-C or T-Hg-T) is held together by
a metal ion rather than by Watson-Crick hydrogen bonds, so the standard
pair recipes, which pull the two bases into hydrogen-bonding geometry, are
wrong for it. This module finds such pairs from the model so that the
restraint generator can leave them unrestrained and the guesser can leave
them out of its canonical guesses.

Rule
----
For each metal ion (a ``HETATM`` atom of a metal element, in a residue that is
not a nucleotide), find the nucleotides that have a metal-binding base atom
(N, O, S or Se, excluding sugar and phosphate atoms) within ``max_distance``
of the ion. A pair is two of those nucleotides that are not neighbours on the
same chain (neighbours are stacked, not paired). If several pairs qualify, the
one with the shortest summed metal-atom distances is chosen and the others
are reported.

A nucleotide is any residue with a C1' atom (or the old-style C1*).

Defaults
--------
Metals are those with atomic number 21 (Sc) or above, as for anomalous
scatterers, so Na, Mg, K and Ca are not considered unless named. The default
distance, 3.6 A, is deliberately generous: refined M-N/O distances are about
2.0-2.4 A, but ions placed by hand in a starting model are often further
from (or clash with) the bases they bridge.

The module is pure Python and reads fixed-column PDB records of the first
model, like :mod:`restraints.anomalous`.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from .anomalous import (
    ATOMIC_NUMBERS,
    DEFAULT_MIN_ATOMIC_NUMBER,
    WATER_RESIDUES,
    element_symbol,
    infer_element,
    parse_element_list,
)


#: Nonmetals and metalloids; every other element is treated as a metal.
NON_METALS = frozenset({
    "H", "D", "HE", "B", "C", "N", "O", "F", "NE", "SI", "P", "S", "CL", "AR",
    "GE", "AS", "SE", "BR", "KR", "SB", "TE", "I", "XE", "AT", "RN",
})
METALS = frozenset(symbol for symbol in ATOMIC_NUMBERS if symbol not in NON_METALS)

#: Default metals for pair detection: metals with Z >= 21.
DEFAULT_PAIR_METALS = frozenset(
    symbol for symbol in METALS if ATOMIC_NUMBERS[symbol] >= DEFAULT_MIN_ATOMIC_NUMBER
)

#: Maximum metal-to-base-atom distance (Angstrom) for a contact.
DEFAULT_MAX_DISTANCE = 3.6

#: Contacts shorter than this are reported as clashes.
CLASH_DISTANCE = 1.8

#: Elements a metal ion can bind to on a base.
LIGAND_ELEMENTS = frozenset({"N", "O", "S", "SE"})

#: Sugar and phosphate atoms, which never count as base atoms.
BACKBONE_ATOMS = frozenset(
    name.replace("'", star)
    for name in (
        "P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P",
        "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
    )
    for star in ("'", "*")
)

ResidueKey = tuple[str, str]  # (chain character as in the PDB, residue number)


@dataclass(frozen=True)
class Atom:
    record: str
    name: str
    element: str | None
    residue_name: str
    chain: str          # raw chain character (" " when blank), as Biopython uses
    resseq: str         # residue number without insertion code, e.g. "12"
    altloc: str
    xyz: tuple[float, float, float]
    line_number: int

    @property
    def residue(self) -> ResidueKey:
        return (self.chain, self.resseq)


def residue_label(key: ResidueKey, residue_name: str = "") -> str:
    chain, resseq = key
    text = f"{chain.strip() or '_'}:{resseq}"
    return f"{text} {residue_name}" if residue_name else text


@dataclass(frozen=True)
class MetalContact:
    residue: ResidueKey
    residue_name: str
    atom_name: str
    distance: float

    @property
    def label(self) -> str:
        return f"{residue_label(self.residue, self.residue_name)} {self.atom_name} {self.distance:.2f} A"


@dataclass(frozen=True)
class MetalMediatedPair:
    metal: Atom
    first: MetalContact
    second: MetalContact

    @property
    def residues(self) -> tuple[ResidueKey, ResidueKey]:
        return (self.first.residue, self.second.residue)

    @property
    def metal_label(self) -> str:
        alt = f" alt {self.metal.altloc}" if self.metal.altloc else ""
        return f"{self.metal.element} {residue_label(self.metal.residue)}{alt}"

    def describe(self) -> str:
        return (
            f"{residue_label(self.first.residue, self.first.residue_name)} <-> "
            f"{residue_label(self.second.residue, self.second.residue_name)} "
            f"via {self.metal_label} "
            f"({self.first.atom_name} {self.first.distance:.2f} A, "
            f"{self.second.atom_name} {self.second.distance:.2f} A)"
        )


def read_atoms(pdb_filename: str | Path) -> list[Atom]:
    """Read ATOM/HETATM records of the first model."""
    atoms: list[Atom] = []
    seen_model = False
    with Path(pdb_filename).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            record = line[0:6]
            if record == "MODEL ":
                if seen_model:
                    break
                seen_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in ("ATOM  ", "HETATM"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                resseq = str(int(line[22:26]))
            except ValueError:
                continue
            residue_name = line[17:20].strip()
            name_field = line[12:16]
            element = element_symbol(line[76:78]) if len(line) > 76 else None
            if element is None:
                element, _ = infer_element(name_field, residue_name)
            atoms.append(
                Atom(
                    record=record.strip(),
                    name=name_field.strip(),
                    element=element,
                    residue_name=residue_name,
                    chain=line[21:22] or " ",
                    resseq=resseq,
                    altloc=line[16:17].strip(),
                    xyz=xyz,
                    line_number=line_number,
                )
            )
    return atoms


def _adjacent(first: ResidueKey, second: ResidueKey) -> bool:
    return first[0] == second[0] and abs(int(first[1]) - int(second[1])) <= 1


class _Grid:
    """Spatial hash so each metal only checks nearby atoms."""

    def __init__(self, atoms: Iterable[Atom], cell: float):
        self.cell = cell
        self.cells: dict[tuple[int, int, int], list[Atom]] = {}
        for atom in atoms:
            self.cells.setdefault(self._key(atom.xyz), []).append(atom)

    def _key(self, xyz):
        return tuple(int(math.floor(value / self.cell)) for value in xyz)

    def near(self, xyz):
        cx, cy, cz = self._key(xyz)
        for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
            yield from self.cells.get((cx + dx, cy + dy, cz + dz), ())


def find_metal_mediated_pairs(
    pdb_filename: str | Path,
    *,
    max_distance: float = DEFAULT_MAX_DISTANCE,
    metals: Iterable[str] | None = None,
) -> tuple[list[MetalMediatedPair], list[str]]:
    """Return ``(pairs, messages)`` for metal-mediated base pairs in the model.

    ``metals`` limits the metal elements (default: :data:`DEFAULT_PAIR_METALS`).
    Messages explain ions that bind fewer than two bases, alternative pairs
    that were not chosen, and contacts short enough to be clashes.
    """
    if max_distance <= 0:
        raise ValueError("max_distance must be positive")
    wanted = parse_element_list(metals) if metals is not None else DEFAULT_PAIR_METALS

    atoms = read_atoms(pdb_filename)

    by_residue: dict[ResidueKey, list[Atom]] = {}
    for atom in atoms:
        by_residue.setdefault(atom.residue, []).append(atom)
    nucleotides = {
        key for key, members in by_residue.items()
        if any(atom.name in ("C1'", "C1*") for atom in members)
    }

    base_atoms = [
        atom for atom in atoms
        if atom.residue in nucleotides
        and atom.name not in BACKBONE_ATOMS
        and atom.element in LIGAND_ELEMENTS
    ]
    metal_atoms = [
        atom for atom in atoms
        if atom.record == "HETATM"
        and atom.element in wanted
        and atom.residue not in nucleotides
        and atom.residue_name.upper() not in WATER_RESIDUES
    ]

    grid = _Grid(base_atoms, max_distance)
    pairs: list[MetalMediatedPair] = []
    messages: list[str] = []
    used_pairs: set[frozenset[ResidueKey]] = set()

    for metal in metal_atoms:
        closest: dict[ResidueKey, MetalContact] = {}
        for atom in grid.near(metal.xyz):
            distance = math.dist(metal.xyz, atom.xyz)
            if distance > max_distance:
                continue
            current = closest.get(atom.residue)
            if current is None or distance < current.distance:
                closest[atom.residue] = MetalContact(
                    residue=atom.residue,
                    residue_name=atom.residue_name,
                    atom_name=atom.name,
                    distance=distance,
                )

        alt = f" alt {metal.altloc}" if metal.altloc else ""
        label = f"{metal.element} {residue_label(metal.residue)}{alt}"
        contacts = sorted(closest.values(), key=lambda contact: contact.distance)

        for contact in contacts:
            if contact.distance < CLASH_DISTANCE:
                messages.append(
                    f"{label} is only {contact.distance:.2f} A from "
                    f"{residue_label(contact.residue, contact.residue_name)} "
                    f"{contact.atom_name}; the ion position clashes and should be refined."
                )

        options = sorted(
            (
                (first, second)
                for first, second in itertools.combinations(contacts, 2)
                if not _adjacent(first.residue, second.residue)
            ),
            key=lambda option: option[0].distance + option[1].distance,
        )

        if not options:
            if contacts:
                bound = ", ".join(contact.label for contact in contacts)
                messages.append(
                    f"{label} binds base atoms of {bound}, but no two non-neighbouring "
                    f"nucleotides; no metal-mediated pair."
                )
            continue

        first, second = options[0]
        key = frozenset((first.residue, second.residue))
        for other_first, other_second in options[1:]:
            messages.append(
                f"{label}: also within {max_distance:g} A of "
                f"{residue_label(other_first.residue, other_first.residue_name)} and "
                f"{residue_label(other_second.residue, other_second.residue_name)}; "
                f"kept the closer pair."
            )
        if key in used_pairs:
            messages.append(f"{label} bridges an already detected pair; not repeated.")
            continue
        used_pairs.add(key)
        pairs.append(MetalMediatedPair(metal=metal, first=first, second=second))

    return pairs, messages


def metal_bound_residues(
    pairs: Sequence[MetalMediatedPair],
) -> dict[ResidueKey, MetalMediatedPair]:
    """Map each residue in a metal-mediated pair to its (first) pair."""
    bound: dict[ResidueKey, MetalMediatedPair] = {}
    for pair in pairs:
        for residue in pair.residues:
            bound.setdefault(residue, pair)
    return bound
