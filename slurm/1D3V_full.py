from clusterrun import createsbatch
import os

# Solver and simulation parameters
Moments = [4, 6, 8]
closure = "ExtGram"
Knudsen = [0.2]#!, 1.0]
sources = ["relaxation_source"]#!, "relaxation_source"]

angle_name = "MinCond"

# SLURM job parameters
threads = 1
nnodes = 1
time = '03:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
    for M in Moments:
        command = f"julia examples/Shock_1D3V_full.jl {M} {closure} {Kn} {source} {angle_name}"
        
        os.system(command)
        # createsbatch(
        #     command, 
        #     nproc=threads, nnodes=nnodes, 
        #     time=time, mem=memory_request, 
        #     output_file=f"out/1D3V_full/1D3V_angles/slurm_output_1D3V_M{M}_closure{closure}_Kn{Kn}_source{source}_anglename{angle_name}.out"
        # )
