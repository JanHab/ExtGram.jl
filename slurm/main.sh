#!/usr/bin/zsh

### Job Parameters
#SBATCH --job-name=main               # Job name
#SBATCH --output=slurm/out/main_out.txt   # Standard output file
#SBATCH --error=slurm/out/main_out.txt             # Standard error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --time=00:05:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --cpus-per-task=4           # Number of CPU cores per task

#Load necessary modules (if needed)
module load Julia

#Your job commands go here
#For example:
julia examples/main.jl --threads=4

#Optionally, you can include cleanup commands here (e.g., after the job finishes)
#For example:
#rm some_temp_file.txt