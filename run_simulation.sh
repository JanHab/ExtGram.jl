#!/bin/bash
#
#SBATCH --job-name=out/VlasovPoisson/TwoStreamInstability/slurm_output_closureExtGram_T25.0_M14_p3_level8.out
#SBATCH --output=out/VlasovPoisson/TwoStreamInstability/slurm_output_closureExtGram_T25.0_M14_p3_level8.out
#
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#
###SBATCH --account=thes1498

if [ -r /usr/local_host/etc/bashrc ]; then
    . /usr/local_host/etc/bashrc
fi

export PATH=$PATH:/home/$USER/bin

module load Julia

julia examples/TwoStreamInstability_clusterrun.jl 14 ExtGram 25.0 8 3 --threads=1
