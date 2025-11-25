#!/usr/bin/zsh
#
### Job Parameters
#SBATCH --job-name=twostream             # Job name
#SBATCH --output=slurm/out/twostream_out.txt   # Standard output file
#SBATCH --error=slurm/out/twostream_out.txt             # Standard error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --time=01:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#
#Load necessary modules (if needed)
module load Julia
#
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
        srun julia examples/TwoStreamInstability_clusterrun.jl --threads=$threads --M=$M --closure=$closure --T_end=$T_end --base_tree_level=$base_tree_level &
    done

#Optionally, you can include cleanup commands here (e.g., after the job finishes)
#For example:
#rm some_temp_file.txt
