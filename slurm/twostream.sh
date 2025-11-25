#!/usr/bin/zsh
<<<<<<< HEAD
#
=======

>>>>>>> 5f0bb95 (clusterrun for twostreaminstability)
### Job Parameters
#SBATCH --job-name=twostream             # Job name
#SBATCH --output=slurm/out/twostream_out.txt   # Standard output file
#SBATCH --error=slurm/out/twostream_out.txt             # Standard error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --time=01:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
<<<<<<< HEAD
#
#Load necessary modules (if needed)
module load Julia
#
=======

#Load necessary modules (if needed)
# module load Julia

>>>>>>> 5f0bb95 (clusterrun for twostreaminstability)
#Your job commands go here
#For example:
threads=4
# export JULIA_NUM_THREADS=$threads
M_vector=(3 5) #!(4 6 8 10 12)
closure="Gram"
T_end=25.0
base_tree_level=8
for M in "${M_vector[@]}"
    do
        # echo "Starting job for M=$M"
        # srun 
<<<<<<< HEAD
        srun julia examples/TwoStreamInstability_clusterrun.jl --threads=$threads --M=$M --closure=$closure --T_end=$T_end --base_tree_level=$base_tree_level &
=======
        julia examples/TwoStreamInstability_clusterrun.jl --threads=$threads --M=$M --closure=$closure --T_end=$T_end --base_tree_level=$base_tree_level &
>>>>>>> 5f0bb95 (clusterrun for twostreaminstability)
    done

#Optionally, you can include cleanup commands here (e.g., after the job finishes)
#For example:
<<<<<<< HEAD
#rm some_temp_file.txt
=======
#rm some_temp_file.txt
>>>>>>> 5f0bb95 (clusterrun for twostreaminstability)
