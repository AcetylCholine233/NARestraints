"""Heavy-atom detection and phenix.refine anomalous-scatterer groups.

This module finds heavy atoms written as ``HETATM`` records in a PDB file and
turns them into ``refinement.refine.anomalous_scatterers`` groups, so that
phenix.refine can refine their f' and f'' against anomalous data.

It is deliberately independent of Biopython and of the residue workbook: it
reads the fixed-column PDB records directly, so its behaviour is easy to check
and it works for metals and other ions that are not nucleotides.

What counts as a heavy atom
---------------------------
By default, any ``HETATM`` atom (waters excluded) whose element has an atomic
number of at least :data:`DEFAULT_MIN_ATOMIC_NUMBER` (21, scandium). That
includes transition metals such as Ag, Hg, Au, Cu, Zn, Ni, Pd and Pt, and
heavier main-group elements such as Se, Br and I. Light ions whose anomalous
signal is usually negligible at typical wavelengths (Na, Mg, K, Ca) and the
light atoms of ligands (C, N, O, P, S, Cl) are not included unless the caller
names them explicitly.

Element identity
----------------
The element column (PDB columns 77-78) is used when it is filled in.
Otherwise the element is inferred: first from an ion residue whose residue and
atom names are both the element symbol (for example ``AG`` / ``AG``), then from
the atom name following the PDB alignment convention, then from the leading
letters of the atom name. Each site records how its element was determined.

Starting values
---------------
Groups are written with ``f_prime = 0`` and ``f_double_prime = 0`` and both
values refined, matching Phenix's own defaults. One group is written per
element, so all sites of one element share one f'/f'' pair.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


_ELEMENT_SYMBOLS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm",
)

#: Upper-case element symbol -> atomic number. Deuterium is treated as H.
ATOMIC_NUMBERS: dict[str, int] = {
    symbol.upper(): number for number, symbol in enumerate(_ELEMENT_SYMBOLS, start=1)
}
ATOMIC_NUMBERS["D"] = 1

#: Default threshold: scandium (Z = 21) and heavier.
DEFAULT_MIN_ATOMIC_NUMBER = 21

#: Residue names treated as water and never scanned.
WATER_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"})

#: The phenix.refine default strategy (Phenix 2.2.1 refinement parameters).
#: ``group_anomalous`` is not part of the default, so it must be added for the
#: anomalous groups to be refined at all.
PHENIX_DEFAULT_STRATEGY = (
    "individual_sites",
    "individual_sites_real_space",
    "individual_adp",
    "occupancies",
)

ELEMENT_COLUMN = "element column"
INFERRED_FROM_ION_RESIDUE = "inferred from ion residue name"
INFERRED_FROM_ATOM_NAME = "inferred from atom name"


@dataclass(frozen=True)
class HeavyAtomSite:
    """One heavy atom found in a HETATM record."""

    element: str          # upper-case symbol, e.g. "AG"
    atom_name: str        # stripped, e.g. "AG"
    residue_name: str     # stripped, e.g. "AG"
    chain: str            # stripped; empty when the chain column is blank
    resseq: str           # stripped residue number, e.g. "101"
    icode: str            # insertion code or ""
    altloc: str           # alternate-location indicator or ""
    occupancy: float | None
    line_number: int      # 1-based line number in the PDB file
    element_source: str   # ELEMENT_COLUMN or one of the INFERRED_* values

    @property
    def resid(self) -> str:
        return f"{self.resseq}{self.icode}"

    @property
    def label(self) -> str:
        chain = self.chain or "_"
        alt = f" alt {self.altloc}" if self.altloc else ""
        return f"{chain}:{self.resid} {self.residue_name} {self.atom_name}{alt}"

    def selection(self) -> str:
        """Phenix atom selection for this atom (all alternate locations)."""
        parts = []
        if self.chain:
            parts.append(f"chain {self.chain}")
        parts.append(f"resid {self.resid}")
        parts.append(f"name {self.atom_name}")
        return " and ".join(parts)


def element_symbol(value: str) -> str | None:
    """Return the upper-case element symbol for ``value``, or None if unknown."""
    text = "".join(ch for ch in value.strip() if ch.isalpha()).upper()
    return text if text in ATOMIC_NUMBERS else None


def parse_element_list(text: str | Iterable[str]) -> frozenset[str]:
    """Parse ``"AG,HG"`` (or an iterable of symbols) into upper-case symbols.

    Raises ValueError for anything that is not an element symbol.
    """
    if isinstance(text, str):
        items = [item for item in text.replace(";", ",").replace(" ", ",").split(",")]
    else:
        items = list(text)

    symbols: set[str] = set()
    for item in items:
        item = item.strip()
        if not item:
            continue
        symbol = item.upper()
        if symbol not in ATOMIC_NUMBERS:
            raise ValueError(f"Not an element symbol: {item!r}")
        symbols.add(symbol)

    if not symbols:
        raise ValueError("No element symbols were given.")
    return frozenset(symbols)


def infer_element(name_field: str, residue_name: str) -> tuple[str | None, str]:
    """Infer an element from the 4-character atom-name field (columns 13-16).

    Returns ``(symbol or None, how it was inferred)``.
    """
    stripped = name_field.strip()
    letters = "".join(ch for ch in stripped if ch.isalpha()).upper()

    # Ion residues such as "AG" / "AG" or "ZN" / "ZN": residue name and atom
    # name agree and are an element symbol. This also resolves "CA" calcium
    # correctly, which the column convention below would read as carbon.
    if letters and letters == residue_name.strip().upper() and letters in ATOMIC_NUMBERS:
        return letters, INFERRED_FROM_ION_RESIDUE

    # PDB convention: a one-letter element is written in column 14 (so column
    # 13 is blank or a digit); a two-letter element fills columns 13-14.
    padded = name_field.ljust(4)
    if padded[0] == " " or padded[0].isdigit():
        candidate = padded[1].upper()
    else:
        candidate = padded[0:2].upper()
    if candidate in ATOMIC_NUMBERS:
        return candidate, INFERRED_FROM_ATOM_NAME

    # Fall back to the leading letters of the stripped name.
    for length in (2, 1):
        candidate = letters[:length]
        if len(candidate) == length and candidate in ATOMIC_NUMBERS:
            return candidate, INFERRED_FROM_ATOM_NAME

    return None, INFERRED_FROM_ATOM_NAME


def _parse_occupancy(line: str) -> float | None:
    try:
        return float(line[54:60])
    except ValueError:
        return None


def find_heavy_atoms(
    pdb_filename: str | Path,
    *,
    elements: Iterable[str] | None = None,
    min_atomic_number: int = DEFAULT_MIN_ATOMIC_NUMBER,
) -> tuple[list[HeavyAtomSite], list[str]]:
    """Find heavy atoms in the HETATM records of the first model.

    If ``elements`` is given, exactly those elements are reported, whatever
    their atomic number. Otherwise every element with atomic number of at
    least ``min_atomic_number`` is reported.

    Returns ``(sites, warnings)``. Warnings report HETATM atoms whose element
    could not be determined; they are skipped rather than guessed.
    """
    wanted = parse_element_list(elements) if elements is not None else None

    sites: list[HeavyAtomSite] = []
    warnings: list[str] = []
    seen_model = False

    path = Path(pdb_filename)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            record = line[0:6]

            if record == "MODEL ":
                if seen_model:
                    break  # only the first model, as for restraint generation
                seen_model = True
                continue
            if record == "ENDMDL":
                break
            if record != "HETATM":
                continue

            residue_name = line[17:20].strip()
            if residue_name.upper() in WATER_RESIDUES:
                continue

            name_field = line[12:16]
            element = element_symbol(line[76:78]) if len(line) > 76 else None
            source = ELEMENT_COLUMN
            if element is None:
                element, source = infer_element(name_field, residue_name)

            atom_name = name_field.strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()

            if element is None:
                warnings.append(
                    f"Line {line_number}: could not determine the element of "
                    f"HETATM {chain or '_'}:{resseq} {residue_name} {atom_name}; skipped."
                )
                continue

            if wanted is not None:
                if element not in wanted:
                    continue
            elif ATOMIC_NUMBERS[element] < min_atomic_number:
                continue

            sites.append(
                HeavyAtomSite(
                    element=element,
                    atom_name=atom_name,
                    residue_name=residue_name,
                    chain=chain,
                    resseq=resseq,
                    icode=line[26:27].strip(),
                    altloc=line[16:17].strip(),
                    occupancy=_parse_occupancy(line),
                    line_number=line_number,
                    element_source=source,
                )
            )

    return sites, warnings


def group_by_element(sites: Sequence[HeavyAtomSite]) -> dict[str, list[HeavyAtomSite]]:
    """Group sites by element, keeping the order in which elements first appear."""
    groups: dict[str, list[HeavyAtomSite]] = {}
    for site in sites:
        groups.setdefault(site.element, []).append(site)
    return groups


def _group_selection(sites: Sequence[HeavyAtomSite]) -> str:
    selections: list[str] = []
    for site in sites:
        selection = site.selection()
        if selection not in selections:  # alternate locations share one selection
            selections.append(selection)
    if len(selections) == 1:
        return selections[0]
    return " or ".join(f"({selection})" for selection in selections)


def anomalous_phil_block(
    sites: Sequence[HeavyAtomSite],
    *,
    include_strategy: bool = True,
    strategy: Sequence[str] = PHENIX_DEFAULT_STRATEGY,
) -> str:
    """Return a ``refine { ... }`` scope for inside ``refinement { ... }``.

    With no sites, return an empty string: a strategy containing
    ``group_anomalous`` without any group makes phenix.refine stop.
    """
    if not sites:
        return ""

    lines = ["  refine {"]

    if include_strategy:
        selected = [item for item in strategy if item != "group_anomalous"]
        selected.append("group_anomalous")
        lines.append(
            "    # NARestraints --anomalous: Phenix default strategy plus group_anomalous."
        )
        lines.append("    strategy = " + " ".join(f"*{item}" for item in selected))

    lines.append("    anomalous_scatterers {")
    for element, members in group_by_element(sites).items():
        escaped = _group_selection(members).replace('"', '\\"')
        lines.extend([
            f"      # {element}: {len(members)} HETATM record(s)",
            "      group {",
            f'        selection = "{escaped}"',
            "        f_prime = 0",
            "        f_double_prime = 0",
            "        refine = *f_prime *f_double_prime",
            "      }",
        ])
    lines.append("    }")
    lines.append("  }")
    return "\n".join(lines)


def describe_sites(sites: Sequence[HeavyAtomSite]) -> list[str]:
    """Human-readable summary lines for the terminal."""
    if not sites:
        return []
    lines = []
    for element, members in group_by_element(sites).items():
        lines.append(f"Anomalous scatterers: {element} x{len(members)}")
        for site in members:
            note = "" if site.element_source == ELEMENT_COLUMN else f" ({site.element_source})"
            lines.append(f"  {site.label}{note}")
    return lines
