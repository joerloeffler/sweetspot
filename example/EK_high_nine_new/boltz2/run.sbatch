#!/bin/bash
#SBATCH --job-name=boltz2
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --partition=low
#SBATCH --gres=shard:30
#SBATCH --nodelist=ai,pika

set -euo pipefail

JOBDIR="$SLURM_SUBMIT_DIR"
cd "$JOBDIR"

eval "$(/z/linux/bin/micromamba shell hook -s posix)"
micromamba activate "/z/bio/biotools/miniconda/envs/boltz2/"

echo "Starting boltz2 on $(hostname) at $(date)"
snakemake --snakefile "$JOBDIR/Snakefile" --directory "$JOBDIR" --cores 16
echo "Finished boltz2 at $(date)"
