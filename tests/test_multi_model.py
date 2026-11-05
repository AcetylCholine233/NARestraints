"""Multi-model PDB files: every step must describe model 1 only."""

import re
from pathlib import Path

import pytest


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def two_model_1bna(path: Path) -> Path:
    """1BNA as model 1; model 2 is a copy renumbered A 201-212, B 113-124.

    Like a file whose two models are different constructs, model 2 has
    residues that do not exist in model 1.
    """
    coords = [
        line for line in (EXAMPLES / "1BNA.pdb").read_text().splitlines()
        if line.startswith("ATOM")
    ]
    shift = {"A": 200, "B": 100}
    model2 = [
        f"{line[:22]}{int(line[22:26]) + shift[line[21]]:>4}{line[26:]}"
        for line in coords
    ]
    lines = ["MODEL        1", *coords, "ENDMDL", "MODEL        2", *model2, "ENDMDL", "END"]
    path.write_text("\n".join(lines) + "\n")
    return path


def _needs_bio():
    pytest.importorskip("Bio")
    pytest.importorskip("openpyxl")
    pytest.importorskip("numpy")


def test_guesser_only_pairs_residues_of_model_one(tmp_path):
    _needs_bio()
    from restraints.guesser import guess_pairs

    pdb = two_model_1bna(tmp_path / "two.pdb")
    candidates, _ = guess_pairs(pdb)
    assert len(candidates) >= 8  # the 12-bp duplex is still found
    for candidate in candidates:
        for residue in (candidate.first, candidate.second):
            assert int(residue.resid) <= 24, f"{residue.chain}:{residue.resid} is from model 2"


def test_builder_stacks_model_one_only(tmp_path):
    _needs_bio()
    from restraints.builder import build_phil_from_pdb

    pdb = two_model_1bna(tmp_path / "two.pdb")
    out = tmp_path / "two.phil"
    build_phil_from_pdb(pdb, [], out)
    stacks = re.findall(
        r"base1 = chain \S+ and resid (\S+)\s+base2 = chain \S+ and resid (\S+)",
        out.read_text(),
    )
    assert len(stacks) == 22
    assert all(int(a) <= 24 and int(b) <= 24 for a, b in stacks)
