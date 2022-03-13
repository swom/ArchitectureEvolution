#!/bin/bash
# Script to run the simulation
#
# Usage, locally:
#
#   ./run_switch_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_switch_loop.sh
#
# Peregrine directives:
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=arc_evo_sim
#SBATCH --output=arc_evo_%j.log

echo "seed: "$1
echo "architecture: "$2
echo "change frequency: " $3
echo "number of total generations: "$4
echo "mutation type: "$5
echo "selection type: "${6}
echo "selection frequency: "${7}

./src/arc_evo -S $1 -N $2 -C $3 -G $4 -m $5 -s ${6} -f ${7}
