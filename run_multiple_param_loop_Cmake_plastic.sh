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
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=3GB
#SBATCH --job-name=run_arc_evo_loop
#SBATCH --output=run_arc_evo_loop.log


# to get arc_evo.pro command interface 
#run arc_evo.exe with command line argument -h 
module load git
module load CMake
module load binutils
mkdir build_test
cd build_test
cmake -DCMAKE_BUILD_TYPE=Release ..
make

declare -a architectures=("1,2,2,2,1")
declare -a max_architectures=("1,2,2,2,1")
declare -a change_freq_As=(0.1)
declare -a change_freq_Bs=(0.01 0.001)
declare -a sel_strs=(0.1 0.5 1 2)
declare -a gen=(500000)
declare -a mut_types=("weights")
declare -a sel_types=("constant")
declare -a change_freq_types=("regular")
declare -a change_sym_types=("asymmetrical")
declare -a adaptation_periods=("on" "off")
declare -a type_of_responses=("plastic" "constitutive")
declare -a sel_freqs=(1)
declare record_top_ind_freq=1000
declare n_observations_reaction_norm=100
declare n_trials=10

for seed in $(seq 1 10)
do
	for arc in "${architectures[@]}"
	do
		for max_arc in "${max_architectures[@]}"
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
								for adaptation_period in "${adaptation_periods[@]}"
								do
									for change_freq_type in "${change_freq_types[@]}"
									do
										for change_sym_type in "${change_sym_types[@]}"
										do
											for change_freq_B in "${change_freq_Bs[@]}"
											do
												for sel_str in "${sel_strs[@]}"
												do
													for type_of_response in "${type_of_responses[@]}"
													do
													echo $seed $arc $max_arc $change_freq_A $gen $mut_type $sel_type $sel_freq $adaptation_period $change_freq_type $change_sym_type $change_freq_B $sel_str $type_of_response
															sbatch ../run_loop_plastic.sh $seed $arc $max_arc $change_freq_A $gen $mut_type $sel_type $sel_freq $record_top_ind_freq $n_observations_reaction_norm $n_trials $adaptation_period $change_freq_type $change_sym_type $change_freq_B $sel_str $type_of_response
													done
												done
											done
										done
									done
								done
							done
						done
					done
				done 
			done
		done
	done
done
