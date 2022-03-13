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


# to get arc_evo.pro command interface 
#run arc_evo.exe with command line argument -h 
module load CMake
module load binutils
mkdir build
cd build
cmake ..
cmake --build . 


declare -a architectures=("1,2,2,2,1")
declare -a change_freq_As=(0)
declare -a gen=(1000000)
declare -a mut_types=("weights")
declare -a sel_types=("sporadic")
declare -a sel_freqs=(100 1000 10000)

for seed in $(seq 1 10)
do
	for arc in "${architectures[@]}"
	do
		for change_freq_A in "${change_freq_As[@]}"
		do
			for gen in "${gen[@]}"
			do	
				for mut_type in "${mut_types[@]}"
				do
					for sel_type in "${sel_types[@]}"
					do
						for sel_freq in "${sel_freqs[@]}"
						do
								echo $seed $arc $change_freq_A $gen $mut_type $sel_type $sel_freq
              
								sbatch ../run_loop_sporadic.sh $seed $arc $max_arc $change_freq_A $gen $mut_type $sel_type $sel_freq
						done
					done
				done
			done 
		done
	done
done
