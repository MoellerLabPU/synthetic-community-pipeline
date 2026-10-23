#!/bin/bash
# Run the pipeline as a SLURM job. This job only runs Snakemake, which then
# submits one job per step (see profiles/slurm/config.yaml), so it needs
# almost no resources but must stay alive until the last step is done.
#
#   sbatch run_pipeline.sh <config.yaml> <working directory> [snakemake options]
#
# Results (results_dir in the config, if relative) and the job logs
# (logs/slurm/) are written inside the working directory.
# Snakemake's own output: snakemake-<jobid>.out in the folder sbatch was run from.
#
#   sbatch run_pipeline.sh config/config.yaml .                 # everything, in the repo
#   sbatch run_pipeline.sh config/config.yaml /workdir/me/run1  # in another folder
#   sbatch run_pipeline.sh config/config.yaml . -n              # dry run
#SBATCH --job-name=snakemake
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=snakemake-%j.out

set -euo pipefail
config=$1
workdir=$2
shift 2                           # anything else goes to snakemake

cd "$SLURM_SUBMIT_DIR"            # the repository, where sbatch was run
eval "$(conda shell.bash hook)"   # make "conda activate" work in a batch job
conda activate synthcom

snakemake \
    --snakefile workflow/Snakefile \
    --profile profiles/slurm \
    --configfile "$config" \
    --directory "$workdir" \
    "$@"
