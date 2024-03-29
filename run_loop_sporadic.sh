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
#SBATCH --cpus-per-task=5	
#SBATCH --mem=1GB
#SBATCH --job-name=arc_evo_sim
#SBATCH --output=arc_evo_%j.log

echo "seed: "$1
echo "architecture: "$2
echo "max_architecture: "$3
echo "change frequency: " $4
echo "number of total generations: "$5
echo "mutation type: "$6
echo "selection type: "${7}
echo "selection frequency: "${8}
echo "record top ind frequency: "$9
echo "number of observations in reaction norms: "${10}
echo "number of trials: "${11}
echo "adaptation period: "${12}
echo "change frequency type: "${13}
echo "selection strength: "${14}
echo "number of mutations: "${15}

./src/arc_evo -S $1 -N $2 -X $3 -C $4 -G $5 -m $6 -s ${7} -f ${8} -R $9 -p ${10} -n ${11} -H ${12} -z ${13} -T ${14} -u${15}
