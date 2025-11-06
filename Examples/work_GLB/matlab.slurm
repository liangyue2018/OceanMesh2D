#!/bin/bash
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=10

# load matlab
module load matlab-2024a-null

# check input
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch matlab.slurm \"<matlab_command>\"" >&2
    exit 1
fi

CMD="$1"

# run script
$PROGLIST -batch "$CMD"
