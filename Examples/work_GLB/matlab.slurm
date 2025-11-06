#!/bin/bash
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=5

# load matlab
module load matlab-2024a-none

# check input
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch matlab.slurm \"<matlab_command>\"" >&2
    exit 1
fi

CMD="$1"

# run script
$PROGLIST -batch "$CMD"
