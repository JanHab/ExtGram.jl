from clusterrun import createsbatch
import os

# Solver and simulation parameters
Moments = [4, 6, 8]
closure = "ExtGram"
Knudsen = [0.1, 1.0, 10.0]
sources = ["relaxation_source", "relaxation_source", "zero_source"]

angle_names = ["Max", "Arc", "Det"]

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
    for M in Moments:
        for angle_name in angle_names:
            command = f"julia examples/1D3VSlab.jl {M} {closure} {Kn} {source} {angle_name}"
            # os.system(command)
            createsbatch(
                command, 
                nproc=threads, nnodes=nnodes, 
                time=time, mem=memory_request, 
                output_file=f"out/1D3D/1D3D_angles/slurm_output_1D3D_M{M}_closure{closure}_Kn{Kn}_source{source}_anglepair{angle_name}.out"
            )
