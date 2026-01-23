#!/bin/bash
#
#SBATCH --job-name=jn
#SBATCH --output=outfile
#
#SBATCH --cpus-per-task=cppt
#SBATCH --ntasks=1
#SBATCH --nodes=numno
#SBATCH --time=aot
#SBATCH --mem=mem_req
#
#SBATCH --account=thes2188

if [ -r /usr/local_host/etc/bashrc ]; then
    . /usr/local_host/etc/bashrc
fi

export PATH=$PATH:/home/$USER/bin

module load Julia

cmd
