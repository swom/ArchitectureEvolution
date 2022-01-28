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
module load Qt5
export CC=g++
export CXX=g++
make clean
qmake arc_evo.pro
make 

#declare -a architectures=("1,2,2,2,1" "1,5,1" "1,5,1,1")
declare -a architectures=("1,2,2,2,1")
declare -a max_architectures=("1,8,8,8,1")
declare -a change_freq=(0 0.005)
declare -a gen=(50000)
declare -a mut_type=("NRduplication" "NRaddition" "addition" "weights_and_activation")

for seed in $(seq 1 10)
do
	for arc in "${architectures[@]}"
	do
		for max_arc in "${max_architectures[@]}"
		do
			for change_freq in "${change_freq[@]}"
			do
				for gen in "${gen[@]}"
				do	
					for mut_type in "${mut_type[@]}"
					do
						  echo $seed
						  echo $arc
						  echo $max_arc
						  echo $change_freq
						  echo $gen
						  echo $mut_type
						  sbatch run_loop.sh $seed $arc $max_arc $change_freq $gen $mut_type
					done
				done
			done
		done
	done
done 
