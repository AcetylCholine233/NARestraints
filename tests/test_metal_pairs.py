import math
import random
import re
from pathlib import Path

import pytest

from restraints.metal_pairs import (
    find_metal_mediated_pairs,
    metal_bound_residues,
    read_atoms,
)


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def record(serial, name, resname, chain, resseq, xyz, element, *, het=False, altloc=" "):
    kind = "HETATM" if het else "ATOM"
    x, y, z = xyz
    return (
        f"{kind:<6}{serial:>5} {name:<4}{altloc}{resname:>3} {chain}{resseq:>4} "
        f"   {x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}"
        f"          {element:>2}"
    )


class Model:
    """Tiny PDB builder: nucleotides are C1' plus chosen base/backbone atoms."""

    def __init__(self):
        self.lines = []

    def nucleotide(self, chain, resseq, atoms, resname="DC"):
        # C1' sits 5 A further out so it never counts as a contact.
        first = next(iter(atoms.values()))
        c1 = tuple(value * 3 for value in first) if any(first) else (9.0, 9.0, 9.0)
        self.lines.append(record(len(self.lines) + 1, "C1'", resname, chain, resseq, c1, "C"))
        for name, xyz in atoms.items():
            element = "P" if name == "P" else name[0]
            self.lines.append(record(len(self.lines) + 1, name, resname, chain, resseq, xyz, element))
        return self

    def metal(self, element, chain, resseq, xyz=(0.0, 0.0, 0.0), resname=None, altloc=" "):
        resname = resname or element
        self.lines.append(record(
            len(self.lines) + 1, element, resname, chain, resseq, xyz, element,
            het=True, altloc=altloc,
        ))
        return self

    def write(self, tmp_path, name="model.pdb"):
        path = tmp_path / name
        path.write_text("\n".join(self.lines) + "\nEND\n")
        return path


def keys(pair):
    return {(chain.strip(), resid) for chain, resid in pair.residues}


def test_c_ag_c_pair_is_found(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N3": (-2.1, 0.0, 0.0)})
        .nucleotide("B", 20, {"N3": (2.1, 0.0, 0.0)})
        .metal("AG", "Z", 1)
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert len(pairs) == 1
    assert keys(pairs[0]) == {("A", "5"), ("B", "20")}
    assert pairs[0].first.atom_name == "N3"
    assert pairs[0].first.distance == pytest.approx(2.1)
    assert "via AG Z:1" in pairs[0].describe()
    assert messages == []


def test_neighbours_on_one_chain_are_not_a_pair(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N7": (-2.1, 0.0, 0.0)}, resname="DG")
        .nucleotide("A", 6, {"N7": (2.1, 0.0, 0.0)}, resname="DG")
        .metal("AG", "Z", 1)
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert pairs == []
    assert "no two non-neighbouring nucleotides" in messages[0]


def test_ion_on_one_base_is_not_a_pair(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N7": (-2.1, 0.0, 0.0)}, resname="DG")
        .nucleotide("B", 20, {"N3": (6.0, 0.0, 0.0)})
        .metal("ZN", "Z", 1)
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert pairs == []
    assert len(messages) == 1


def test_phosphate_sugar_and_carbon_atoms_do_not_count(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"OP1": (-2.0, 0.0, 0.0), "O4'": (0.0, -2.0, 0.0)})
        .nucleotide("B", 20, {"OP2": (2.0, 0.0, 0.0), "C5": (0.0, 2.0, 0.0)})
        .metal("AG", "Z", 1)
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert pairs == [] and messages == []


def test_light_metals_ignored_by_default_and_found_when_named(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N7": (-2.1, 0.0, 0.0)}, resname="DG")
        .nucleotide("B", 20, {"O6": (2.1, 0.0, 0.0)}, resname="DG")
        .metal("MG", "Z", 1)
        .write(tmp_path)
    )
    assert find_metal_mediated_pairs(pdb)[0] == []
    pairs, _ = find_metal_mediated_pairs(pdb, metals="MG")
    assert keys(pairs[0]) == {("A", "5"), ("B", "20")}


def test_waters_and_atom_records_are_not_metals(tmp_path):
    lines = (
        Model()
        .nucleotide("A", 5, {"N3": (-2.1, 0.0, 0.0)})
        .nucleotide("B", 20, {"N3": (2.1, 0.0, 0.0)})
        .lines
    )
    lines.append(record(90, "AG", "AG", "Z", 1, (0, 0, 0), "AG"))  # ATOM, not HETATM
    pdb = tmp_path / "m.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n")
    assert find_metal_mediated_pairs(pdb)[0] == []


def test_closest_pair_is_chosen_and_alternatives_reported(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N3": (-2.1, 0.0, 0.0)})
        .nucleotide("B", 20, {"N3": (2.2, 0.0, 0.0)})
        .nucleotide("C", 40, {"O2": (0.0, 3.3, 0.0)})
        .metal("AG", "Z", 1)
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert keys(pairs[0]) == {("A", "5"), ("B", "20")}
    assert any("kept the closer pair" in message for message in messages)


def test_distance_cutoff_and_clash_report(tmp_path):
    pdb = (
        Model()
        .nucleotide("B", 4, {"O4": (1.19, 0.0, 0.0)}, resname="DU")
        .nucleotide("A", 12, {"N4": (-3.52, 0.0, 0.0)})
        .metal("AG", "B", 101)
        .write(tmp_path)
    )
    # The generous default catches a poorly placed starting-model ion...
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert keys(pairs[0]) == {("B", "4"), ("A", "12")}
    assert any("clashes" in message for message in messages)
    # ...while a strict cutoff does not.
    assert find_metal_mediated_pairs(pdb, max_distance=3.0)[0] == []
    with pytest.raises(ValueError):
        find_metal_mediated_pairs(pdb, max_distance=0)


def test_alternate_locations_do_not_duplicate_a_pair(tmp_path):
    pdb = (
        Model()
        .nucleotide("A", 5, {"N3": (-2.1, 0.0, 0.0)})
        .nucleotide("B", 20, {"N3": (2.1, 0.0, 0.0)})
        .metal("AG", "Z", 1, (0.0, 0.1, 0.0), altloc="A")
        .metal("AG", "Z", 1, (0.0, -0.1, 0.0), altloc="B")
        .write(tmp_path)
    )
    pairs, messages = find_metal_mediated_pairs(pdb)
    assert len(pairs) == 1
    assert any("already detected" in message for message in messages)
    bound = metal_bound_residues(pairs)
    assert set(bound) == {("A", "5"), ("B", "20")}


def test_grid_search_matches_brute_force(tmp_path):
    rng = random.Random(7)
    model = Model()
    for i in range(60):
        xyz = tuple(rng.uniform(-12, 12) for _ in range(3))
        model.nucleotide("ABCD"[i % 4], 10 * i, {"N3": xyz})
    metals = [tuple(rng.uniform(-12, 12) for _ in range(3)) for _ in range(40)]
    for j, xyz in enumerate(metals):
        model.metal("HG", "Z", j + 1, xyz)
    pdb = model.write(tmp_path)

    atoms = read_atoms(pdb)
    bases = [a for a in atoms if a.name == "N3"]
    pairs, _ = find_metal_mediated_pairs(pdb, max_distance=4.0)
    found = {frozenset(keys(pair)) for pair in pairs}

    # Brute force: for each metal, the closest pair of bases within 4 A.
    # (No residues are chain neighbours: residue numbers step by 10.)
    expected = set()
    for xyz in metals:
        near = sorted(
            (math.dist(a.xyz, xyz), (a.chain.strip(), a.resseq))
            for a in bases if math.dist(a.xyz, xyz) <= 4.0
        )
        if len(near) >= 2:
            expected.add(frozenset((near[0][1], near[1][1])))
    assert found == expected
    assert len(found) >= 5  # the random model really exercises the search


def test_bundled_dfpucker_silver_bridges_df_and_dg():
    pairs, _ = find_metal_mediated_pairs(EXAMPLES / "dfPucker.pdb")
    assert [keys(pair) for pair in pairs] == [{("A", "9"), ("D", "8")}]
    pair = pairs[0]
    assert (pair.first.residue_name, pair.first.atom_name) == ("DF", "S1")
    assert (pair.second.residue_name, pair.second.atom_name) == ("DG", "N2")


# ---------------------------------------------------------------------------
# Builder and guesser behaviour (need Biopython and openpyxl).


def _needs_bio():
    pytest.importorskip("Bio")
    pytest.importorskip("openpyxl")


def _stacking_pairs(text):
    return re.findall(
        r"base1 = chain (\S+) and resid (\S+)\s+base2 = chain (\S+) and resid (\S+)", text
    )


def test_waters_and_ions_are_never_stacked(tmp_path):
    _needs_bio()
    from restraints.builder import build_phil_from_pdb

    out = tmp_path / "1bna.phil"
    build_phil_from_pdb(EXAMPLES / "1BNA.pdb", [], out)
    stacks = _stacking_pairs(out.read_text())
    # 1BNA: two 12-nt strands -> 11 stacks each; its 80 waters share chain A.
    assert len(stacks) == 22
    assert all(int(r1) <= 24 and int(r2) <= 24 for _, r1, _, r2 in stacks)


def test_builder_skips_all_restraints_on_metal_mediated_bases(tmp_path):
    _needs_bio()
    from restraints.base_pairs import read_base_pair_file
    from restraints.builder import build_phil_from_pdb

    pairs = tmp_path / "pairs.txt"
    pairs.write_text("A 9\nD 8\n")
    out = tmp_path / "df.phil"
    build_phil_from_pdb(EXAMPLES / "dfPucker-fixed.pdb", read_base_pair_file(pairs), out)
    text = out.read_text()

    assert not re.search(r"chain A and resid 9\b", text)
    assert not re.search(r"chain D and resid 8\b", text)
    # No new stack bridges the removed base.
    assert ("A", "8", "A", "10") not in _stacking_pairs(text)

    kept = tmp_path / "df-keep.phil"
    build_phil_from_pdb(
        EXAMPLES / "dfPucker-fixed.pdb", read_base_pair_file(pairs), kept,
        keep_metal_pair_stacking=True,
    )
    kept_text = kept.read_text()
    assert re.search(r"chain A and resid 9\b", kept_text)  # stacking only
    assert "distance_ideal" not in kept_text  # still no base-pair restraints

    old = tmp_path / "df-old.phil"
    build_phil_from_pdb(EXAMPLES / "dfPucker-fixed.pdb", [], old, metal_pairs=False)
    assert re.search(r"chain A and resid 9\b", old.read_text())


def test_guesser_never_guesses_a_metal_mediated_pair(tmp_path):
    _needs_bio()
    from restraints.base_pairs import read_base_pair_file
    from restraints.guesser import guess_pairs_detailed, write_guess

    candidates, _, metal_pairs = guess_pairs_detailed(
        EXAMPLES / "dfPucker-fixed.pdb", allow_noncanonical=True
    )
    assert [keys(pair) for pair in metal_pairs] == [{("A", "9"), ("D", "8")}]
    for candidate in candidates:
        for residue in (candidate.first, candidate.second):
            assert (residue.chain, residue.resid) not in {("A", "9"), ("D", "8")}

    out = tmp_path / "guessed.txt"
    write_guess(out, candidates, metal_pairs)
    text = out.read_text()
    assert "# Metal-mediated pair" in text
    # The commented block is ignored when the file is used as a pair file.
    assert len(read_base_pair_file(out)) == len(candidates)
