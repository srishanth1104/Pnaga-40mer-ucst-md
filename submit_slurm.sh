#!/usr/bin/env bash
#SBATCH --job-name=pnaga40_ucst
#SBATCH --output=logs/pnaga_%A_%a.out
#SBATCH --error=logs/pnaga_%A_%a.err
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
# Uncomment and set your cluster's GPU partition name below:
##SBATCH --partition=gpu

set -euo pipefail
mkdir -p logs

# Uncomment and adjust to match your cluster's GROMACS GPU module:
# module load gromacs/2024-gpu

TEMPS=(30 40 50 60 70 80)
TEMP_C=${TEMPS[$SLURM_ARRAY_TASK_ID]}

bash 02_run_ucst_gpu.sh "${TEMP_C}"
