import os
import subprocess
import sys
from pathlib import Path

import pytest

from restraints.readyset import (
    READYSET_OPTIONS,
    ReadySetError,
    find_ready_set,
    prepare_model,
)


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

FAKE_PROGRAM = r'''
import os, sys
from pathlib import Path
mode = os.environ.get("FAKE_READYSET", "copy")
model = Path(sys.argv[1])
Path("fake_args.txt").write_text("\n".join(sys.argv[1:]) + "\n")
print("fake phenix.ready_set", " ".join(sys.argv[1:]))
if mode == "fail":
    print("Sorry: something went wrong")
    sys.exit(1)
if mode == "noop":
    print("No unknown residues")
    sys.exit(0)
if mode == "silent":
    sys.exit(0)
lines = model.read_text().splitlines()
if mode == "hydrogen":
    lines.insert(-1, "HETATM 9999  H1  HOH W 999       0.000   0.000   0.000  1.00 20.00           H")
if mode == "drop":
    lines = [l for l in lines if not l.startswith("HETATM")] + ["END"]
Path(model.stem + ".updated.pdb").write_text("REMARK   fake ready_set\n" + "\n".join(lines) + "\n")
Path(model.stem + ".ligands.cif").write_text("data_comp_list\n")
'''


def make_fake_ready_set(folder: Path) -> Path:
    """A phenix.ready_set stand-in. A shell wrapper avoids shebang paths with spaces."""
    folder.mkdir(parents=True, exist_ok=True)
    script = folder / "fake_ready_set.py"
    script.write_text(FAKE_PROGRAM)
    wrapper = folder / "phenix.ready_set"
    wrapper.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n')
    wrapper.chmod(0o755)
    return wrapper


def hetatm(serial, name, resname, chain, resseq, element, xyz=(0.0, 0.0, 0.0)):
    x, y, z = xyz
    return (
        f"HETATM{serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00          {element:>2}"
    )


def write_model(folder: Path, name="model.pdb", *, metal=True) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1"]
    lines.append(hetatm(1, "O", "HOH", "W", 1, "O", (5.0, 0.0, 0.0)))
    if metal:
        lines.append(hetatm(2, "AG", "AG", "Z", 6, "AG"))
    lines.append("END")
    path = folder / name
    path.write_text("\n".join(lines) + "\n")
    return path


class fake_mode:
    """Set the stand-in's behaviour for one block."""

    def __init__(self, mode):
        self.mode = mode

    def __enter__(self):
        self.old = os.environ.get("FAKE_READYSET")
        os.environ["FAKE_READYSET"] = self.mode

    def __exit__(self, *exc):
        if self.old is None:
            os.environ.pop("FAKE_READYSET", None)
        else:
            os.environ["FAKE_READYSET"] = self.old


def test_model_without_heavy_atoms_is_left_alone(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    model = write_model(tmp_path / "data", metal=False)
    result = prepare_model(model, executable=fake)
    assert not result.ran and result.model == model
    assert not (tmp_path / "data" / "model_readyset").exists()


def test_heavy_atom_triggers_readyset_with_required_options(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    data = tmp_path / "学习一个" / "Structural Biology"  # spaces and non-ASCII
    model = write_model(data)
    with fake_mode("copy"):
        result = prepare_model(model, executable=fake)

    workdir = data / "model_readyset"
    assert result.ran and "AG x1" in result.reason
    assert result.directory == workdir.resolve()
    assert result.model == workdir.resolve() / "model.updated.pdb"
    assert result.model.read_text().startswith("REMARK   fake ready_set")
    assert result.ligand_cif == workdir.resolve() / "model.ligands.cif"
    assert result.log.is_file() and not result.noop and result.warnings == []

    args = (workdir / "fake_args.txt").read_text().split()
    assert args[0] == "model.pdb"  # a copy, by relative name, inside workdir
    for option in READYSET_OPTIONS:
        assert option in args
    assert "ready_set.actions.hydrogens=False" in args
    assert not any("hydrogens=True" in arg for arg in args)
    assert "ready_set.actions.ligands=True" in args

    # Nothing but the work folder is written beside the input model.
    assert sorted(p.name for p in data.iterdir()) == ["model.pdb", "model_readyset"]


def test_always_and_never_modes(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    plain = write_model(tmp_path / "a", metal=False)
    with fake_mode("copy"):
        assert prepare_model(plain, mode="always", executable=fake).ran
    heavy = write_model(tmp_path / "b")
    result = prepare_model(heavy, mode="never", executable=tmp_path / "missing")
    assert not result.ran and result.model == heavy
    with pytest.raises(ValueError):
        prepare_model(heavy, mode="sometimes")


def test_noop_uses_input_and_ignores_stale_outputs(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    model = write_model(tmp_path / "data")
    stale_dir = tmp_path / "data" / "model_readyset"
    stale_dir.mkdir()
    (stale_dir / "model.updated.pdb").write_text("stale\n")
    with fake_mode("noop"):
        result = prepare_model(model, executable=fake)
    assert result.ran and result.noop
    assert result.model == model.resolve()
    assert not (stale_dir / "model.updated.pdb").exists()


def test_added_hydrogens_are_refused(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    model = write_model(tmp_path / "data")
    with fake_mode("hydrogen"), pytest.raises(ReadySetError, match="hydrogen"):
        prepare_model(model, executable=fake)


def test_failure_and_missing_output_are_errors(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    model = write_model(tmp_path / "data")
    with fake_mode("fail"), pytest.raises(ReadySetError, match="exit status 1"):
        prepare_model(model, executable=fake)
    with fake_mode("silent"), pytest.raises(ReadySetError, match="wrote no"):
        prepare_model(model, executable=fake)


def test_changed_atom_count_is_warned(tmp_path):
    fake = make_fake_ready_set(tmp_path / "phenix")
    model = write_model(tmp_path / "data")
    with fake_mode("drop"):
        result = prepare_model(model, executable=fake)
    assert any("changed the atom count" in warning for warning in result.warnings)


def test_find_ready_set_from_phenix_folder_and_environment(tmp_path):
    install = tmp_path / "phenix-2.2.1-6174"
    fake = make_fake_ready_set(install / "bin")
    assert find_ready_set(install) == fake
    assert find_ready_set(fake) == fake
    with pytest.raises(ReadySetError, match="not an executable"):
        find_ready_set(tmp_path / "nothing-here")

    old = os.environ.get("NARESTRAINTS_READY_SET")
    os.environ["NARESTRAINTS_READY_SET"] = str(install)
    try:
        assert find_ready_set() == fake
    finally:
        if old is None:
            os.environ.pop("NARESTRAINTS_READY_SET")
        else:
            os.environ["NARESTRAINTS_READY_SET"] = old


def test_command_line_runs_readyset_before_restraints(tmp_path):
    pytest.importorskip("Bio")
    pytest.importorskip("openpyxl")
    fake = make_fake_ready_set(tmp_path / "phenix")
    data = tmp_path / "data"
    data.mkdir()
    model = data / "1bna-ag.pdb"
    lines = [l for l in (EXAMPLES / "1BNA.pdb").read_text().splitlines() if not l.startswith("END")]
    lines.append(hetatm(9000, "AG", "AG", "Z", 6, "AG", (50.0, 50.0, 50.0)))
    model.write_text("\n".join(lines) + "\nEND\n")
    pairs = data / "pairs.txt"
    pairs.write_text("")

    env = dict(os.environ, FAKE_READYSET="copy")
    completed = subprocess.run(
        [sys.executable, "-m", "restraints", str(model), str(pairs),
         "--anomalous", "--ready-set", str(fake)],
        cwd=str(Path(__file__).resolve().parents[1]),
        env=env, stdin=subprocess.DEVNULL, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    assert completed.returncode == 0, completed.stdout
    assert "ReadySet: run because of heavy atoms" in completed.stdout
    assert "1bna-ag.updated.pdb" in completed.stdout
    assert (data / "1bna-ag_readyset" / "fake_args.txt").is_file()
    phil = (data / "1bna-ag.phil").read_text()
    assert 'selection = "chain Z and resid 6 and name AG"' in phil

    # --readyset never skips it.
    completed = subprocess.run(
        [sys.executable, "-m", "restraints", str(model), str(pairs),
         "--readyset", "never", "-o", str(data / "never.phil")],
        cwd=str(Path(__file__).resolve().parents[1]),
        env=env, stdin=subprocess.DEVNULL, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    assert completed.returncode == 0, completed.stdout
    assert "ReadySet: skipped (--readyset never)" in completed.stdout
