# HELIUS project on oral and gut microbiomes
Interdisciplinary project within the UvA Research Priority Area - Personal Microbiome Health, using HELIUS oral and gut microbiome data

## Data
Two microbiome data types from the HELIUS cohort, processed in parallel:
- **16S rRNA amplicon sequencing** — throat and nose swabs
- **Shotgun metagenomics** — tongue and throat swabs

Raw sequencing output and clinical/metadata are cleaned into `phyloseq`/table objects in `data/processed/` (see `scripts/1a_datacleaning_helius.R` and `scripts/1b_datacleaning_biome.R`).

## Analysis pipeline
Scripts in `scripts/` are numbered in run order:

| Script | Analysis |
|---|---|
| `1a_datacleaning_helius.R` | Clean HELIUS clinical/metadata |
| `1b_datacleaning_biome.R` | Clean and filter 16S and shotgun microbiome data into `phyloseq` objects |
| `2_tableone.R` | Table 1: cohort characteristics |
| `3_relative_abundance_plots.R` | Compositional (stacked bar) plots of taxon relative abundance |
| `4_alpha_diversity_16s_throat.R` | Alpha diversity, 16S throat (descriptive) |
| `5_alpha_diversity_16s_ethnicity.R` | Alpha diversity, 16S throat and nose, stratified by ethnicity |
| `6_beta_diversity_16s.R` | Beta diversity, 16S throat and nose, stratified by ethnicity (PERMANOVA, betadisper) |
| `7_beta_diversity_16s_migration.R` | Beta diversity, 16S, non-Dutch groups pooled by migration generation/acculturation |
| `8_differential_abundance_16s.R` | Differential abundance, 16S throat and nose, pairwise ethnicity comparisons (MaAsLin2) |
| `9_upset_diffabund_16s.R` | Overlap of significant differentially abundant taxa across ethnicity pairs |
| `10_alpha_diversity_shotgun.R` | Alpha diversity, shotgun tongue and throat, stratified by ethnicity |

Shotgun beta diversity and differential abundance analyses, analogous to the 16S ones above, are planned but not yet implemented.

Outputs (plots, tables) are written to `results/`, grouped by analysis type.

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