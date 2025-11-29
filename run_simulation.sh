#!/bin/bash
#
#SBATCH --job-name=out/Riemann1D/slurm_output_Kn1.0_T0.3_N2000_p1_level8.out
#SBATCH --output=out/Riemann1D/slurm_output_Kn1.0_T0.3_N2000_p1_level8.out
#
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#
###SBATCH --account=thes1498

if [ -r /usr/local_host/etc/bashrc ]; then
    . /usr/local_host/etc/bashrc
fi

export PATH=$PATH:/home/$USER/bin

module load Julia

julia examples/Riemann1DBGK_clusterrun.jl 2000 0.3 8 1 1.0
