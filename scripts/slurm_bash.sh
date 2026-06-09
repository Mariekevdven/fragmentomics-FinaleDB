#!/bin/bash

#SBATCH --job-name=frag_extraction
#SBATCH --output=logs/frag_extraction_%A_%a.out
#SBATCH --error=logs/frag_extraction_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000%10

mkdir -p logs output

## source "$(conda info --base)/etc/profile.d/conda.sh"
## conda activate fragmentomics

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/processing_manifest.txt)

python scripts_marieke/Slurm_script.py \
    --input "$INPUT" \
    --reference-genome /home/d.vessies@NKI/projects/reference_genome/oncoref_out/GRCh38_oncoref_v1.fa \
    --outdir output/

