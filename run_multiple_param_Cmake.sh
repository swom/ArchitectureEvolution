#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_param_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_param_loop.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_arc_evo_loop
#SBATCH --output=run_arc_evo_loop.log
module load Qt5
module load GCC
module load CMake
cmake -S . -B build
cmake --build . --config Release --target install