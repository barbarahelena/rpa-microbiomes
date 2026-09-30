# Readability refactor validation

Baseline: `3704c99ef3b9018e88219d9b272f5a86bab6a536`.
Updated: 2026-09-29. The full pipeline is intentionally left for the user to run.

## Current bounded checks

Run from the repository root with Python 3.11+ and the existing Pixi environment installed:

```bash
python3 tests/run_readability_checks.py
```

This runner does not invoke Pixi tasks, access cohort data, or run the full
pipeline. It reads the baseline from Git and creates isolated temporary folders.
It checks:

- R parsing, preserved public task names, task dependencies and helper paths;
  unchanged package dependencies and lockfile.
- Identical parsed cleaning expressions, shared theme, ethnicity beta-diversity
  computation and statistical helper bodies against the baseline.
- Missing-cache errors for all eight cache-backed reports, with an identified
  producer and no attempt to recompute.
- Migration computation/reporting on two synthetic sites, retaining the existing
  999 permutations and eligibility threshold. All 30 baseline CSV/PDF outputs
  match; PDF comparisons ignore only date metadata. The report also runs with
  processed inputs removed and leaves the saved caches unchanged.
- ASV and genus reporting with no significant taxa, one significant taxon and
  multiple significant taxa. These tests deliberately stub MaAsLin2 with fixed
  output tables to check the structural boundary, not inference. Model input
  objects/settings, summary rows, RNG state and generated reports match baseline.
- One qualifying ethnicity group: both differential-abundance scripts skip fitting
  and plots. This also checks the no-pair UpSet reporting path.

The final bounded suite passed. A separate static check found one producer for
11 representative computation caches, exported tables and native model outputs.
Alpha-report output patterns explicitly include only CSV/PDF files, excluding
RDS computation caches. `git diff --check` and Pixi task discovery also passed.

The approved plotting-seed handoff retains the random draw formerly consumed
by each labelled volcano plot between model fits. It does not introduce a new
fixed production seed. Tests cover both the draw and no-draw branches.

## Earlier full-data comparisons (2026-09-28)

The previous session ran full-data cohort, exposure, sample-characteristics,
relative-abundance, alpha-diversity, geography and manuscript-figure comparisons;
reduced-data ethnicity beta-diversity and ASV/genus comparisons; and full ASV
compute/report plus full genus computation comparisons.

Across completed comparisons, 1,163 CSV/TSV files, 313 PDFs and 23 PNGs matched
exactly (ignoring only PDF date metadata). Figure 1 and its supplements matched
in both formats. No differences were found in completed comparisons. These
counts are the previous session's recorded results, not reruns on 2026-09-29.

An isolated Pixi test verified report cache hits, rebuilding after report edits,
unchanged computation-cache content/mtime, preserved compute cache hits after
report edits, and rebuilding after compute edits. It used `--skip-deps` with
prepared inputs, not a fresh full-pipeline execution.

The earlier `/tmp/rpa-readability-baseline` and `/tmp/rpa-readability-validation`
files are no longer present. The checked-in tests above reconstruct their own
baseline from Git and print the new temporary artifact/log directory.

Earlier comparisons used the installed environment, external seed `20260928`,
one beta-diversity worker, and `DIFF_AB_TEST_N=40` / `BETA_DIV_TEST_N=40` for reduced
runs. Geography checks required `PROJ_DATA` to point to the installed environment's
`share/proj` directory in both versions; no package or analysis fix was made.

## Full-run validation left to the user

Fresh cleaning, full-data ethnicity beta-diversity rebuilding, full migration
comparison and completion of full-genus reporting validation remain outstanding.
The synthetic checks do not establish full-cohort result equivalence.

The refactor changes script paths and adds compute/report handoffs, so the first
normal task run may rebuild prerequisites and create missing caches. Run the
pipeline with test-mode overrides unset when ready:

```bash
pixi run pipeline
```

Historical unseeded permutation p-values are not guaranteed to reproduce exactly;
their existing production RNG behaviour is unchanged. Keep any historical outputs
needed for comparison before rerunning, since normal runs reuse their filenames.

Existing repository data/results were not overwritten by these checks. Changes
remain uncommitted.
