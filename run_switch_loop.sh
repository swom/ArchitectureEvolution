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
#SBATCH --job-name=single_switch_sim
#SBATCH --output=switch_%j.log

echo "seed: "$1
echo "architecture: "$2
./mutational_switches --net_arc $2 --seed $1 
