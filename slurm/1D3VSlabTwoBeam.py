from clusterrun import createsbatch
import os

# Solver and simulation parameters
Moments = [4, 6, 8]
closure = "ExtGram"
Kn = 1.0
source = "relaxation_source"
angle_name = "Arc"

# SLURM job parameters
threads = 1
nnodes = 1
time = '01:00:00'
memory_request = '16G'

# Two Beam
for M in Moments:
    command = f"julia --project=. examples/TwoBeam_1D3VSlab.jl {M} {closure} {Kn} {source} {angle_name}"
    os.system(command)
    # createsbatch(
    #     command, 
    #     nproc=threads, nnodes=nnodes, 
    #     time=time, mem=memory_request, 
    #     output_file=f"out/1D3V/1D3V_TwoBeam/slurm_output_1D3V_M{M}_closure{closure}_Kn{Kn}_source{source}_anglepair{angle_name}.out"
    # )
