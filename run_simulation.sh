#!/bin/bash
#
#SBATCH --job-name=out/VlasovPoisson/LandauDamping/slurm_output_closureExtGram_T25.0_M25_p3_level6.out
#SBATCH --output=out/VlasovPoisson/LandauDamping/slurm_output_closureExtGram_T25.0_M25_p3_level6.out
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

julia examples/LandauDamping_clusterrun.jl 25 ExtGram 25.0 6 3
