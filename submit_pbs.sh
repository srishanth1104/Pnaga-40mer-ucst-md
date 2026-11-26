#!/usr/bin/env bash
#PBS -N pnaga40_ucst
#PBS -l select=1:ncpus=8:ngpus=1:mem=16gb
#PBS -l walltime=48:00:00
#PBS -j oe
# Uncomment and set your PBS GPU queue name below:
##PBS -q gpu

set -euo pipefail
cd "${PBS_O_WORKDIR}"

# Uncomment and adjust to match your cluster's GROMACS GPU module:
# module load gromacs

MAX_HOURS=47.5 bash 02_run_ucst_gpu.sh
