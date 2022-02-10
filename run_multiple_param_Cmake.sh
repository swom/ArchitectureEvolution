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

module load CMake
module load binutils
mkdir build
cd build
cmake ..
cmake --build . 

declare -a architectures=("2,2,2,2,1")
declare -a max_architectures=("2,8,8,8,1")
declare -a change_freq_A=(0.005)
declare -a gen=(1000000)
declare -a mut_types=("NRaddition" "NRduplication")
declare -a change_types=("regular")

for seed in $(seq 1 10)
do
	for arc in "${architectures[@]}"
	do
		for max_arc in "${max_architectures[@]}"
		do
			for change_freq_A in "${change_freq_A[@]}"
			do
				for gen in "${gen[@]}"
				do	
					for mut_type in "${mut_types[@]}"
					do
            for change_type in "${change_types[@]}"
            do
						   echo "seed:  "$seed
						   echo "architecture:  "$arc
						   echo "max_architecture:  "$max_arc
						   echo "change_freq_A:  "$change_freq_A
						   echo "number of generations:  "$gen
						   echo "mutation type:  "$mut_type
						   echo "environmental change type:  "$change_type
              
						  sbatch ../run_loop.sh $seed $arc $max_arc $change_freq_A $gen $mut_type $change_type
					done
				done
			done
		done
	done
done 
done