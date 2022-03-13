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
echo "max_arc: " $3
echo "change frequency: " $4
echo "number of total generations: "$5
echo "mutation type: "$6
echo "change type:  "$7
echo "duplication rate: "$8
echo "activation rate: "$9  
echo "selection type: "${10}
echo "selection frequency: "${11}

./src/arc_evo -S $1 -N $2 -X $3 -C $4 -G $5 -m $6 -e$7 -D $8 -A $9 -s ${10} -f ${11}
