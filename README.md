# HELIUS project on oral and gut microbiomes
Interdisciplinary project within the UvA Research Priority Area - Personal Microbiome Health, using HELIUS oral and gut microbiome data

## Data
Two microbiome data types from the HELIUS cohort, processed in parallel:
- **16S rRNA amplicon sequencing** — throat and nose swabs
- **Shotgun metagenomics** — tongue and throat swabs

Raw sequencing output and clinical/metadata are cleaned into `phyloseq`/table objects in `data/processed/` (see `scripts/1a_datacleaning_helius.R` and `scripts/1b_datacleaning_biome.R`).

### Upstream read processing
`scripts/0_run_vsearch.sh` is a SLURM batch script that runs a Nextflow vsearch pipeline (clustering/rarefaction of raw sequencing reads into ASV tables) on the Snellius HPC cluster, ahead of and separate from the `pixi`-managed steps below. Its output feeds into `1a_datacleaning_helius.R`/`1b_datacleaning_biome.R`; submit it with `sbatch scripts/0_run_vsearch.sh`.

## Analysis pipeline
Scripts in `scripts/` are numbered in run order:

| Script | Analysis | Pixi task |
|---|---|---|
| `1a_datacleaning_helius.R` | Clean HELIUS clinical/metadata | `clean-helius` |
| `1b_datacleaning_biome.R` | Clean and filter 16S and shotgun microbiome data into `phyloseq` objects | `clean-biome` |
| `2_tableone.R` | Table 1: cohort characteristics | `tableone` |
| `3_airpollution_participants.R` | Air pollution exposure distribution among HELIUS participants, overall and by ethnicity | `airpollution-participants` |
| `4_airpollution_amsterdam.R` | Amsterdam-wide PC6 air pollution map | `airpollution-amsterdam` |
| `5_relative_abundance_plots.R` | Compositional (stacked bar) plots of taxon relative abundance | `relabund-plots` |
| `6_alpha_diversity_16s_ethnicity.R` | Alpha diversity, 16S throat and nose, stratified by ethnicity | `alpha-16s` |
| `7a_beta_diversity_16s_compute.R` | Beta diversity, 16S throat and nose: permutation-heavy PERMANOVA/betadisper computation, cached to `.rds` | `beta-16s-compute` |
| `7b_beta_diversity_16s_report.R` | Beta diversity, 16S throat and nose: rebuild plots/tables from cached PERMANOVA/betadisper results (PCoA, betadisper, PERMANOVA, covariate screen, ethnicity attenuation) | `beta-16s-report` |
| `8_beta_diversity_16s_migration.R` | Beta diversity, 16S, non-Dutch groups pooled by migration generation/acculturation | `beta-16s-migration` |
| `9_differential_abundance_16s.R` | Differential abundance, 16S throat and nose, pairwise ethnicity comparisons (MaAsLin2) | `diffabund-16s` |
| `9b_differential_abundance_genus_16s.R` | Differential abundance at genus level, 16S throat and nose, pairwise ethnicity comparisons (MaAsLin2) | `diffabund-16s-genus` |
| `10_upset_diffabund_16s.R` | Overlap of significant differentially abundant taxa across ethnicity pairs | `upset-diffabund-16s` |
| `11_alpha_diversity_shotgun.R` | Alpha diversity, shotgun tongue and throat, stratified by ethnicity | `alpha-shotgun` |

Shotgun beta diversity and differential abundance analyses, analogous to the 16S ones above, are planned but not yet implemented.

Each script can be run individually via its pixi task, e.g. `pixi run beta-16s`, or the whole pipeline via `pixi run pipeline`. `pixi task list` shows all tasks with descriptions.

Outputs (plots, tables) are written to `results/`, grouped by analysis type.

### Reusing results and detecting stale caches

Analysis tasks use [Pixi's file-content cache](https://pixi.prefix.dev/v0.66.0/workspace/advanced_tasks/#caching).
Run tasks normally, for example `pixi run figure2` or `pixi run pipeline`.
Pixi checks prerequisites first and prints `cache hit` when an analysis can
be skipped. An existing result is reusable only after a successful cached
task run and while its recorded inputs and outputs still match.

Each task declares its own `inputs` and `outputs` in `pixi.toml`:

- Input **contents** include the analysis script, the files it reads, and
  `pixi.lock`. Editing covariates, filtering rules or thresholds in a script
  therefore invalidates that task. Touching a file without changing its
  contents does not invalidate it.
- Missing or modified output files invalidate their producer. Output
  patterns are scoped to the producer, even when several tasks share a
  results directory. For example, writing beta-diversity reports does not
  invalidate the separate distance/model cache.
- R and the analysis packages are resolved to exact builds in `pixi.lock`, so
  a lockfile change invalidates the analysis tasks. `Maaslin2` is the sole
  exception: Bioconda has no R 4.4 build, so the setup script installs it with
  BiocManager.

For example, changing only `13_figure2_assembly.R` rebuilds Figure 2 while
reusing valid cleaned data and beta-diversity models. Changing raw metadata
reruns cleaning; downstream tasks rerun when their input contents change.
Changing the locked software environment invalidates all analyses.

`beta-16s-report` now validates its compute dependency too. If that cache is
stale, the compute task runs before reporting; `beta-16s` is an alias for
this same dependency chain. Cache validation applies to normal **Pixi task
runs**. Direct `Rscript` calls and `pixi run --skip-deps` bypass parts of
this protection and should not be used to establish that results are current.
Tasks use `Rscript --vanilla` so user `.Rprofile`/`.Renviron` files cannot
silently change an analysis behind the cache's back. Export supported
settings in the shell instead.

The first run after enabling caching recomputes the requested tasks and
their prerequisites: old result files have no trusted cache record. A
failed task does not establish a successful cache record. Pixi 0.66 has no
`--force` option for cached tasks. To rerun one analysis deliberately, run
its command directly in the Pixi environment, for example:

```bash
pixi run Rscript --vanilla scripts/13_figure2_assembly.R
```

This bypasses task dependencies and does not refresh the task's cache record;
use it only when those dependencies are already current. Alternatively,
removing or changing a declared output makes the normal task run stale. To
inspect the exact files Pixi hashes, use `pixi run -vvv figure2` (this executes
stale tasks, not just a status check).

When adding an input file or helper script, declare it in the task's inputs;
undeclared inputs cannot be detected automatically. Environment-controlled
test settings and changes to the separately installed `Maaslin2` package are
also not cache inputs. Use the direct command above when intentionally using
one of those overrides. Caching preserves a successful run; it does not by
itself make stochastic analyses reproducible across rebuilds. Script 7a
explicitly seeds its computations as described below.

### Resuming beta-diversity computation (7a)

Script 7a has three explicit settings at the top:

```r
BASELINE_R2_PERMUTATIONS <- 0L
INFERENCE_PERMUTATIONS <- 999L
BETA_DIV_SEED <- 42L
```

The attenuation baselines contribute only R², which does not require a
permutation test. Setting their permutation count to zero avoids computing
unused p-values. All reported p-values and covariate-selection tests still
use 999 permutations. If baseline p-values become a requested output, raise
the baseline setting and add those p-values to the output explicitly.

Each independent test gets a deterministic seed derived from the master
seed, site, distance metric and job name. Thus skipping a completed test,
changing worker count or changing baseline permutations does not shift the
random sequence used by another test. Newly seeded p-values can differ from
historical unseeded results; repeated runs of the new code are reproducible
within the same software environment.

Internal checkpoints live under `results/beta_diversity/checkpoints/`
(or `results/beta_diversity_test/checkpoints/`). They cover distances,
individual covariate/pair tests, shared complete-case R² baselines,
dispersion fits/tests, completed PERMANOVA blocks and PCoA ordinations.
Every checkpoint validates the site's input-file contents, computation
code, settings and software versions, plus a checksum of its saved result.
A change to only the nose input preserves throat
checkpoints. Worker count alone does not invalidate internal checkpoints.
Changed code conservatively invalidates all 7a checkpoints.

Checkpoints are published only after successful completion using an atomic
file replacement. If interrupted, restart with the same normal command:

```bash
pixi run pipeline
# Or finish just beta-diversity computation:
pixi run beta-16s-compute
```

Look for `Checkpoint hit:` and `Computing:` messages identifying the site,
metric and test. Completed jobs are reused; interrupted jobs are recomputed.
The original per-site `.rds` files keep their reporting fields and are
assembled atomically after both distances finish. Old files from before this
checkpoint implementation cannot establish a resumable run.

## Setup

### Prerequisites
Install [pixi](https://pixi.sh):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Install the environment
Clone the repository and run the setup script:

```bash
git clone <repo-url>
cd rpa-microbiomes
chmod +x setup_pixi.sh
./setup_pixi.sh
```

This installs all dependencies (phyloseq, tidyverse, vegan, decontam, etc.) and additional Bioconductor packages.

### Run R in the pixi environment
```bash
pixi run R
```

Or activate the environment shell:
```bash
pixi shell
R
```

## Contributors
Roel van der Ploeg, Kevin Singh, Barbara Verhaar

## Funding
This project was supported by an RPA-PMH seed grant
