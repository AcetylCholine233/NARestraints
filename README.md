# NARestraints

The [reviewed workbook audit](docs/workbook-audit-2026-10-11.md) documents
unreleased 5BU/BRU metadata, DF/S6G atom-name corrections and the authorized
canonical-C5 whitespace cleanup. Tagged v1.1.3 retains its original workbook.

For maintainers, document substantial features for human readers with purpose,
inputs/outputs, method, tools, evidence and limitations, including a methods
paragraph when useful for publication. Update that feature description whenever
its behavior or dependencies change; keep planned capabilities explicit.

Generate Phenix nucleic-acid restraint PHIL from a PDB model and a reviewed
base-pair range file. NARestraints also provides a base-pair guesser and a D/L
nucleic-acid mirroring tool. It is usable independently of NASolve.

## Install version 1.1.3

Use a dedicated Python environment (Python 3.10 or later):

```bash
python3 -m venv .venv
.venv/bin/python -m pip install "git+https://github.com/vecchioni-lab/NARestraints.git@v1.1.3"
```

The GitHub release also supplies a wheel and source archive. Install the wheel
with `python -m pip install /path/to/narestraints-1.1.3-py3-none-any.whl`
inside the intended environment. This is a GitHub release, not a PyPI upload.
The `restraints` import namespace and bundled `restraints/data/Ligands.xlsx`
remain unchanged. Phenix and Coot are separate applications, not bundled here.

For development, use `python -m pip install -e '.[test]'` from the checkout.
When NASolve uses an editable NARestraints checkout, install it with NASolve's
own interpreter, and check imports from inside NASolve or outside both projects.
Running an import check inside NARestraints can hide an older installed copy.

## Command-line tools

From an activated environment, using explicit input paths:

```bash
narestraints input.pdb pairs.txt -o restraints.phil
narestraints-guesser input.pdb -o guessed-pairs.txt
narestraints-mirror input.pdb -o mirrored.pdb
```

`python -m restraints` is also supported. Run each command with `--help` for
options. The guesser accepts `--non-canonical` to opt into its existing
noncanonical recipes. Review guessed pairs before using them for refinement.

A range file contains two chain/range lines per paired stretch:

```text
A 103:108
C 214:209

A 109:115
B 125:119
```

Blank lines separate stretches; ranges may run forward or backward.

## Anomalous scatterers (`--anomalous`)

**Status: unreleased feature branch.** Add `--anomalous` to also write the
model's heavy atoms as phenix.refine anomalous-scatterer groups, so their f'
and f'' are refined against anomalous data:

```bash
narestraints input.pdb pairs.txt -o restraints.phil --anomalous
narestraints input.pdb pairs.txt -o restraints.phil --anomalous-elements AG
```

**Inputs and selection.** Only `HETATM` records of the first model are read;
waters (`HOH`, `WAT`, `DOD` and similar) are skipped. By default an atom is a
heavy atom when its element has atomic number 21 (Sc) or more: this covers
transition metals such as Ag, Hg, Au, Cu, Zn, Ni, Pd and Pt, and Se, Br and I,
including those inside modified nucleotides written as `HETATM`. Na, Mg, K, Ca
and ligand C/N/O/P/S/Cl are not included by default; name them with
`--anomalous-elements` (for example `--anomalous-elements AG,MG`), which
replaces the cutoff with an explicit list.

**Element identity.** The element column (PDB columns 77-78) is used when
present. If it is blank, the element is inferred from an ion residue whose
residue and atom names are the same element symbol (`AG`/`AG`), then from the
PDB atom-name alignment convention, then from the atom name's leading letters.
Inferred elements are reported in the terminal; undeterminable ones are
skipped with a warning.

**Output.** One group per element, inside `refinement.refine`:

```text
  refine {
    strategy = *individual_sites *individual_sites_real_space *individual_adp *occupancies *group_anomalous
    anomalous_scatterers {
      group {
        selection = "(chain Z and resid 6 and name AG) or (chain Z and resid 7 and name AG)"
        f_prime = 0
        f_double_prime = 0
        refine = *f_prime *f_double_prime
      }
    }
  }
```

f' and f'' start at 0 (Phenix's defaults) and are both refined. Alternate
locations of one atom share one selection. `group_anomalous` is not part of
Phenix's default strategy, so the file sets the default strategy plus
`group_anomalous`. Because a PHIL strategy line replaces the whole strategy,
use `--no-anomalous-strategy` when the strategy is set elsewhere (another
PHIL file, the command line, or NASolve) and add `group_anomalous` there.
Without `--anomalous`, the output is byte-for-byte unchanged.

**Evidence.** Pure-Python tests cover element detection, inference, cutoffs,
alternate locations, insertion codes, blank chains, multi-model files and the
bundled `dfPucker.pdb` (Ag), `2j02.pdb` (Zn, Mg) and `1BNA.pdb` (water only)
examples. The generated scope was checked against the Phenix 2.2.1
`phenix/refinement/__init__.params` definitions with Phenix's own PHIL parser:
every parameter was recognised, and a misspelled scope was rejected. A native
`phenix.refine` run has not yet been recorded.

**Limitations.** Refining f'' requires anomalous observations (F+/F- or
I+/I-) to be supplied to phenix.refine. Phenix does not allow anomalous group
refinement with experimental phases under the MLHL target. Starting values are
not calculated from the wavelength. `ATOM` records are never scanned. Only the
restraint generator gains this option; the guesser and mirror tools are
unchanged.

## Metal-mediated base pairs and stacking

**Status: unreleased feature branch.** Both `narestraints` and
`narestraints-guesser` look for metal-mediated base pairs (for example
C-Ag-C or T-Hg-T) in the model and treat them differently from hydrogen-bonded
pairs.

**Detection.** A metal ion is a `HETATM` atom of a metal element (default:
atomic number 21/Sc and above; `--metal-elements AG,HG` to choose) in a
residue that is not a nucleotide and not water. A nucleotide is any residue
with a C1' atom. For each ion, NARestraints finds nucleotides with a base
N, O, S or Se atom (not sugar or phosphate atoms) within `--metal-distance`
(default 3.6 A) of the ion. Two such nucleotides that are not neighbours on
one chain form the metal-mediated pair; if several qualify, the pair with the
shortest summed distances is used and the others are reported. The default
distance is generous because ions placed in starting models are often too far
from, or clashing with, the bases they bridge; refined M-N/O distances are
about 2.0-2.4 A. Contacts under 1.8 A are reported as clashes.

**Restraint generation.** Any pair-file pair that involves a base of a
metal-mediated pair gets no pair restraints (H-bond distances, angles,
parallelity or planarity), with a warning if the listed partner differs from
the detected one. These bases also get no stacking restraints unless
`--keep-metal-pair-stacking` is given; their neighbours are not stacked across
the gap.

**Guessing.** Bases of metal-mediated pairs are never candidates, never extend
a run and never fill a gap. Detected pairs are printed and written to the
guessed pair file as commented blocks, which `narestraints` ignores.

**Stacking of non-nucleotides.** Independently of the above, stacking is only
generated between nucleotide-like residues (workbook-mapped, or with a C1'
atom). Ions, waters and ligands that share a chain with the nucleic acid are
no longer stacked onto the neighbouring base. This changes output for models
whose chains contain such residues.

`--no-metal-pairs` turns detection off and restores the previous pair and
guessing behaviour; the non-nucleotide stacking fix always applies.

**Evidence and limitations.** Pure-Python tests cover detection rules, cutoffs,
clashes, alternate locations and a brute-force check of the spatial search;
builder and guesser tests use `examples/dfPucker.pdb`, where Ag bridges A:9
DF S1 (2.22 A) and D:8 DG N2 (2.55 A), and `examples/1BNA.pdb`. Detection is
purely geometric: it does not check coordination chemistry, and an ion that
bridges two bases without forming a designed pair is treated the same way.
Insertion codes are ignored when matching residues.

## Reviewed recipe library

Pair identity and restraint geometry are separate. Pair lookup is order-independent;
the mapped roles within a recipe are ordered.

| Pair | Recipe | First role / second role |
| --- | --- | --- |
| A:T | `AT` | A / T |
| G:C | `GC` | G / C |
| D:T | `D_T` | D / T |
| B:S | `GC` | B / S |
| Z:P | `GC` | P / Z |
| K:X | `GC` | X / K |
| I:C | `AT` | I / C |
| A:G | `AG_IX` | A / G |
| G:T | `GU_XXVIII` | G / T |

`AG_IX` and `GU_XXVIII` currently provide distance restraints without angle
arrays; the guesser requires `--non-canonical` to consider them. Among recognized
residue categories, a pair with no registered recipe produces a warning and
skips pair-specific restraints. Missing mappings or unsupported categories can
still stop generation; the `Other` extension is planned work. The explicit table
in `restraints/recipe_library.py` is the implementation authority; the bundled
atom mappings are in `restraints/data/Ligands.xlsx`.

## What changed in 1.1.3

B:S remains B=G-like/S=C-like; the reversed **Z:P** and **K:X** recipe roles
are corrected to **P=G-like/Z=C-like** and **X=G-like/K=C-like**, with the
same orientation used by modified-base stacking planes. These pairs retain
their independent identities; the shared GC geometry template is not
a `force = G:C` substitution. The distinct D:T recipe remains unchanged.
The committed `Ligands.xlsx`, atom mappings, monomer chemistry and
previous v1.1.2 release are not altered.

Workbook-tab coverage now verifies B:S, Z:P, K:X, D:T and reverse orientation.
A native NASolve Z:P preparation and subsequent numerically successful
refinement have been demonstrated using an independent NASolve W-template
compatibility correction. This does **not** prove the experimental chemistry
of the donor reflections or mark independent B:S/K:X native campaigns complete.

See [v1.1.3 release notes](docs/release-v1.1.3.md) for release evidence and
limitations. The [current development handoff](docs/development-handoff.md)
records later integration progress and open work: the combined WC-like native
preparation/linkage/coverage diagnostics passed, while geometry qualification
and a fresh full-auto acceptance run remain pending. Historical
[v1.1.2 release notes](docs/release-v1.1.2.md) and the
[archived handoff/port notes](docs/history/README.md) preserve prior evidence.

## Release verification

`bash scripts/verify_release.sh` builds wheel and source distributions, installs
each into a fresh environment, runs the full test suite from a separate directory,
and checks all three console tools and the bundled workbook. It is the Linux
release/CI harness; normal package use is not restricted to Linux.

The release workflow tests Python 3.10, 3.12, and 3.14. A push changing
`pyproject.toml` on `main` can publish **v1.1.3** only after those tests
pass; a manual dispatch on `main` can retry that same gated release. Existing tags and
releases are not overwritten. See the Actions result for execution status.
