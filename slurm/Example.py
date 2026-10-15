from clusterrun import createsbatch
import os

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

command = f"julia examples/FirstVlasovMaxwell.jl"
# os.system(command)
createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out"
)
