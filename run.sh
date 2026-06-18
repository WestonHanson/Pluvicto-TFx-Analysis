#!/bin/bash
#SBATCH --job-name=test_sbatch_run
#SBATCH --output=pon_%A_%a.out
#SBATCH --error=pon_%A_%a.err
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --partition=campus-new

ml Python/3.12.3-GCCcore-13.3.0

python /fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/main.py