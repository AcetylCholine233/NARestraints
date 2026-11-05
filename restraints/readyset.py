"""Prepare a model with phenix.ready_set before restraint generation.

NARestraints and phenix.refine both need chemical definitions for every
residue in the model. Models with heavy atoms (metal ions, iodinated or
brominated nucleotides and similar) usually also contain HETATM residues that
need ligand restraints, so by default such models are first passed through
phenix.ready_set, and the updated model it writes is used instead.

ReadySet is always run with

* ``ready_set.actions.hydrogens=False``: never add hydrogens;
* ``ready_set.actions.ligands=True``: always generate ligand restraints;
* ``ready_set.actions.optimise_ligand_geometry=False``: leave coordinates alone.

It runs in its own directory (default ``<model stem>_readyset`` next to the
input) on a copy of the input, so nothing is written beside the original
model. A result that contains hydrogens is refused. If ReadySet reports that
there was nothing to do and writes no updated model, the input is used
unchanged.

The phenix.ready_set executable is looked for in this order: an explicit path,
the ``NARESTRAINTS_READY_SET`` environment variable, ``PATH``, ``$PHENIX``,
and standard install locations (``~/phenix-*``, ``/Applications/phenix-*``,
``/usr/local/phenix-*``, ``/opt/phenix-*``).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from .anomalous import element_symbol, find_heavy_atoms, infer_element


READYSET_OPTIONS = (
    "ready_set.actions.hydrogens=False",
    "ready_set.actions.ligands=True",
    "ready_set.actions.optimise_ligand_geometry=False",
)

MODES = ("auto", "always", "never")

PROGRAM = "phenix.ready_set"


class ReadySetError(RuntimeError):
    """ReadySet could not be found, failed, or produced an unusable model."""


@dataclass
class ReadySetResult:
    model: Path                 # the model to use from now on
    ran: bool                   # whether phenix.ready_set was run
    reason: str                 # why it was or was not run
    directory: Path | None = None
    log: Path | None = None
    ligand_cif: Path | None = None
    command: tuple[str, ...] = ()
    noop: bool = False          # ReadySet ran but wrote no updated model
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Deciding whether to run


def heavy_atom_reason(pdb_filename: str | Path, elements: Iterable[str] | None = None) -> str | None:
    """Return a short description of the heavy atoms that trigger ReadySet, or None."""
    sites, _ = find_heavy_atoms(pdb_filename, elements=elements)
    if not sites:
        return None
    counts: dict[str, int] = {}
    for site in sites:
        counts[site.element] = counts.get(site.element, 0) + 1
    listed = ", ".join(f"{element} x{count}" for element, count in counts.items())
    return f"heavy atoms in HETATM records ({listed})"


# ---------------------------------------------------------------------------
# Finding phenix.ready_set


def _version_key(path: Path):
    parts = re.findall(r"\d+", path.name)
    return tuple(int(part) for part in parts)


def _candidates_in_install(root: Path) -> list[Path]:
    return [root / "bin" / PROGRAM, root / "phenix_bin" / PROGRAM, root / "build" / "bin" / PROGRAM]


def find_ready_set(explicit: str | Path | None = None) -> Path:
    """Locate phenix.ready_set; raise ReadySetError with advice if not found."""

    def usable(path: Path) -> bool:
        return path.is_file() and os.access(path, os.X_OK)

    if explicit:
        path = Path(explicit).expanduser()
        if path.is_dir():
            for candidate in _candidates_in_install(path):
                if usable(candidate):
                    return candidate
            raise ReadySetError(f"No executable {PROGRAM} found in Phenix folder {path}")
        if usable(path):
            return path
        raise ReadySetError(f"{path} is not an executable {PROGRAM}")

    env = os.environ.get("NARESTRAINTS_READY_SET")
    if env:
        return find_ready_set(env)

    on_path = shutil.which(PROGRAM)
    if on_path:
        return Path(on_path)

    phenix = os.environ.get("PHENIX")
    if phenix:
        for candidate in _candidates_in_install(Path(phenix)):
            if usable(candidate):
                return candidate

    roots: list[Path] = []
    for base in (Path.home(), Path("/Applications"), Path("/usr/local"), Path("/opt")):
        try:
            roots.extend(p for p in base.glob("phenix-*") if p.is_dir())
        except OSError:
            continue
    for root in sorted(roots, key=_version_key, reverse=True):
        for candidate in _candidates_in_install(root):
            if usable(candidate):
                return candidate

    raise ReadySetError(
        f"Could not find {PROGRAM}. Pass --ready-set /path/to/phenix-VERSION "
        "(the Phenix folder) or /path/to/phenix.ready_set, set NARESTRAINTS_READY_SET, "
        "or use --readyset never to skip ReadySet."
    )


# ---------------------------------------------------------------------------
# Running it


def _atom_counts(pdb_filename: Path) -> tuple[int, int, int]:
    """Return (atoms, HETATM atoms, hydrogens) over all coordinate records."""
    atoms = hetatms = hydrogens = 0
    with pdb_filename.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[0:6]
            if record not in ("ATOM  ", "HETATM"):
                continue
            atoms += 1
            hetatms += record == "HETATM"
            element = element_symbol(line[76:78]) if len(line) > 76 else None
            if element is None:
                element, _ = infer_element(line[12:16], line[17:20])
            if element in ("H", "D"):
                hydrogens += 1
    return atoms, hetatms, hydrogens


_STALE_PATTERNS = ("*.updated.pdb", "*.ligands.cif", "ready_set.log")


def run_ready_set(
    pdb_filename: str | Path,
    executable: str | Path,
    directory: str | Path | None = None,
    *,
    runner: Callable = subprocess.run,
) -> ReadySetResult:
    """Run phenix.ready_set on a copy of the model and return the result."""
    source = Path(pdb_filename).resolve()
    if not source.is_file():
        raise ReadySetError(f"Model not found: {source}")

    workdir = Path(directory) if directory else source.with_name(f"{source.stem}_readyset")
    workdir = workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    for pattern in _STALE_PATTERNS:
        for stale in workdir.glob(pattern):
            if stale.resolve() != source:  # never delete the input itself
                stale.unlink()

    # Run on a copy, by relative name, inside the work directory: ReadySet
    # writes next to its working directory, and this keeps long or
    # non-ASCII folder names out of its arguments.
    local_name = source.name if source.suffix.lower() == ".pdb" else f"{source.stem}.pdb"
    local = workdir / local_name
    if local.resolve() != source:
        shutil.copyfile(source, local)

    command = (str(executable), local_name, *READYSET_OPTIONS)
    log = workdir / "ready_set.log"
    try:
        completed = runner(
            list(command),
            cwd=str(workdir),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError as exc:
        raise ReadySetError(f"Could not start {executable}: {exc}") from exc
    output = completed.stdout or ""
    log.write_text(output, encoding="utf-8")

    if completed.returncode != 0:
        tail = "\n".join(output.strip().splitlines()[-8:])
        raise ReadySetError(
            f"{PROGRAM} failed (exit status {completed.returncode}); see {log}\n{tail}"
        )

    stem = Path(local_name).stem
    updated = workdir / f"{stem}.updated.pdb"
    if not updated.is_file():
        others = sorted(workdir.glob("*.updated.pdb"))
        updated = others[0] if len(others) == 1 else updated

    result = ReadySetResult(
        model=updated,
        ran=True,
        reason="",
        directory=workdir,
        log=log,
        command=command,
    )

    if not updated.is_file():
        if re.search(r"no unknown residues", output, re.IGNORECASE) or re.search(
            r"nothing to (do|be done)", output, re.IGNORECASE
        ):
            result.model = source
            result.noop = True
        else:
            raise ReadySetError(
                f"{PROGRAM} finished but wrote no {stem}.updated.pdb; see {log}"
            )

    cifs = sorted(workdir.glob("*.ligands.cif"))
    if cifs:
        expected = workdir / f"{stem}.ligands.cif"
        result.ligand_cif = expected if expected in cifs else cifs[0]

    before = _atom_counts(source)
    after = _atom_counts(result.model)
    if after[2] > before[2]:
        raise ReadySetError(
            f"{PROGRAM} added {after[2] - before[2]} hydrogen atom(s) despite "
            f"hydrogens=False; refusing to use {result.model}"
        )
    if after[:2] != before[:2]:
        result.warnings.append(
            f"{PROGRAM} changed the atom count: {before[0]} atoms ({before[1]} HETATM) "
            f"before, {after[0]} ({after[1]} HETATM) after. Check {result.model}."
        )
    return result


def prepare_model(
    pdb_filename: str | Path,
    *,
    mode: str = "auto",
    executable: str | Path | None = None,
    directory: str | Path | None = None,
    heavy_elements: Iterable[str] | None = None,
    runner: Callable = subprocess.run,
) -> ReadySetResult:
    """Decide whether to run ReadySet, run it if needed, and return the model to use.

    ``mode`` is ``"auto"`` (run when the model has heavy atoms in HETATM
    records), ``"always"`` or ``"never"``.
    """
    if mode not in MODES:
        raise ValueError(f"ReadySet mode must be one of {MODES}, not {mode!r}")
    source = Path(pdb_filename)

    if mode == "never":
        return ReadySetResult(model=source, ran=False, reason="--readyset never")

    if mode == "auto":
        reason = heavy_atom_reason(source, heavy_elements)
        if reason is None:
            return ReadySetResult(
                model=source, ran=False, reason="no heavy atoms in HETATM records"
            )
    else:
        reason = "--readyset always"

    result = run_ready_set(source, find_ready_set(executable), directory, runner=runner)
    result.reason = reason
    return result


def report(result: ReadySetResult, echo: Callable[[str], None] = print) -> None:
    """Print what happened, for the command-line tools."""
    if not result.ran:
        if result.reason.startswith("no heavy"):
            return  # ordinary models: stay quiet
        echo(f"ReadySet: skipped ({result.reason}).")
        return
    echo(f"ReadySet: run because of {result.reason}.")
    echo(f"ReadySet: {' '.join(result.command)}")
    echo(f"ReadySet: log {result.log}")
    for warning in result.warnings:
        echo(f"WARNING: {warning}")
    if result.noop:
        echo("ReadySet: nothing to change; using the input model unchanged.")
    else:
        echo(f"ReadySet: using the updated model {result.model}")
        echo("ReadySet: refine this updated model; the restraints refer to it.")
    if result.ligand_cif is not None:
        echo(
            f"ReadySet: ligand restraints {result.ligand_cif} "
            "(give this file to phenix.refine as well)"
        )
