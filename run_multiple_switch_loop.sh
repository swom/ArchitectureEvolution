#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_switches.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_switches.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_multiple_switches
#SBATCH --output=run_multiple_switches.log


# simulation_logic_only has this command-linae interface:
#
# simulation_logic_only --[function] s[seed] N[net_architecture]
module load Qt5
export CC=g++
export CXX=g++
make clean
qmake mutational_switches.pro
make 

declare -a architectures=("1,2,1" "1,5,1" "1,5,1,1")

for i in $(seq 1 10)
do
  for j in "${architectures[@]}"
do
  echo $i
  echo $j
  sbatch run_switch_loop.sh $i $j
  done
done 
