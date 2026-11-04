# Changelog

## Unreleased

- Detect metal-mediated base pairs from coordinates (`restraints/metal_pairs.py`):
  a HETATM metal (Z >= 21 by default) within 3.6 A of base N/O/S/Se atoms of
  two non-neighbouring nucleotides. `narestraints` skips pair restraints and,
  unless `--keep-metal-pair-stacking`, stacking for those bases;
  `narestraints-guesser` never guesses them and records them as comments in
  its output. Options `--metal-distance`, `--metal-elements`,
  `--no-metal-pairs`. New `guess_pairs_detailed()`; `guess_pairs()` keeps its
  two-value return.
- Stop generating stacking between nucleotides and following non-nucleotide
  residues (ions, waters, ligands) in the same chain.
- Add an opt-in `--anomalous` / `--anomalous-elements` toggle to `narestraints`
  that finds heavy atoms in HETATM records (atomic number >= 21 by default,
  waters excluded) and writes one phenix.refine
  `refinement.refine.anomalous_scatterers` group per element, with f'/f''
  starting at 0 and refined, plus the default strategy with `group_anomalous`
  (suppress with `--no-anomalous-strategy`). New pure-Python module
  `restraints/anomalous.py` and `tests/test_anomalous.py`. Output without the
  toggle is unchanged. PHIL checked against Phenix 2.2.1 parameter definitions;
  no native phenix.refine run recorded yet.
- Correct swapped 5BU/BRU sugar types and abbreviations, including the sorting
  view; align reviewed DF O2 mapping with atom S1, and trim all 129 literal
  `C5 ` role mappings. Document the exact 138 changed cells, source evidence
  and before/after workbook hashes; preserve unrelated values and workbook parts.
  Pair targets, numerical dictionaries and protonation are unchanged.
- Add workbook-backed regressions and document the remaining coupled PSU
  glycosidic/base-plane issue. This is unreleased source data; installed v1.1.3
  remains separate until explicitly updated.

- Refresh the active documentation checkpoint against v1.1.3 and the subsequent
  native NASolve combined-run receipts, retaining the boundary between
  preparation/linkage/coverage passes and pending geometry/full-auto validation.
- Put the planned SMILES/CIF audit after the geometry test and before Topo/Scout
  in the NASolve-owned workflow. No audit implementation or scientific code changes.
- Consolidate current recipe documentation in the README; preserve the superseded
  port README and September handoff as exact, checksummed historical backups.

## 1.1.3 — 2026-10-08

- Preserve the existing named B:S, Z:P, K:X and D:T recipes and their
  experimental identities. NARestraints' `GC` is the shared three-contact
  *geometry template*, not a user-forced G:C identity substitution.
- Correct the reversed Z:P and K:X roles: B/P/X are guanine-like (N1/N2/O6)
  and S/Z/K cytosine-like (O2/N3/N4). Keep B:S unchanged.
- Align modified-base stacking-plane selection with those role assignments.
- Replace misleading Z/K guanine-like synthetic fixtures and add tests using
  actual `Ligands.xlsx` records from the B, S, Z, P, K, X and D sheets,
  plus a thymine partner for the separate `D_T` recipe. Also verify pair order
  reversal and the exact contact/plane atom selections.
- D denotes the reviewed diaminopurine-like role in D:T; no D:A recipe is
  registered, and none is invented in this patch. I:C stays unchanged pending
  a separate evidence-based assessment; all other recipes and the ligand
  workbook are untouched.
- GitHub candidate CI passed Python 3.10/3.12/3.14 (source and installed
  wheel/sdist) before release version bump; final 1.1.3 CI must also pass.
  The user reported 21/21 focused and a green corrected full NARestraints
  regression from the independent candidate worktree.
- The corrected Z:P pair passed native NASolve PostMR with three checked
  P:Z contacts and no geometry-force override. After an independent
  NASolve Saenger-template compatibility correction for modified off-pair
  bases, native Phenix AutoRefine reached numerical SOLVED/refine-001 with
  frozen input integrity OK. User reported good visual bonds and planes in
  Coot (mild plane deviations acceptable).
- Real B:S/K:X native NASolve tests with unavailable local component
  dictionaries remain pending and are **not** represented as passing.
  The synthetic modified-base W challenge is not proof of the true
  experimental chemistry. Existing historical runs and the dirty local
  NARestraints feature checkout remain untouched.

## 1.1.2 — 2026-09-15

- Release the previously merged modified-residue stacking fix (PR #1): explicit
  parallelity restraints for mapped modified neighbours, native stacking for
  canonical neighbours, and a warning-preserving fallback for unknown mappings.
- Bump the distribution version from 1.1.1 so the release is distinguishable
  from the older installation that lacked the stacking module.
- Preserve the dictionary-audit handoff and add an authoritative current release
  checkpoint rather than losing historical reasoning during branch cleanup.
- Add wheel/source-distribution build and non-editable installation checks,
  including all console tools and the unchanged committed ligand workbook.
- No changes to pair recipes, sulfur target values, or the ligand workbook.
  NASolve's 1AP/OP3 integration and unfinished DE repair are separate projects.

## 1.1.1 — 2026-08-25

Installable package, bundled ligand workbook, working-directory-independent
resource lookup, and CLI entry points for restraints, guessing, and mirroring.
