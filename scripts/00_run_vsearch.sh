#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=50:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate nextflow
mkdir -p /projects/prjs1784/rpa/results_260130
nextflow run /projects/prjs1784/pipelines/vsearchpipeline/main.nf \
        -profile snellius \
        --input /projects/prjs1784/rpa/data/sample_sheet.csv \
        --outdir /projects/prjs1784/rpa/results_260130 \
        --skip_primers \
        --cluster_minsize 10 \
        --cluster_alpha 1.5 \
        --rarelevel 10000 \
        --email 'bjh.verhaar@gmail.com' \
        -resume
