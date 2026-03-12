from clusterrun import createsbatch
import os

# Solver and simulation parameters
Moments = [4, 8, 14, 5, 9, 15]
closure = "ExtGram"
T_end = 25.0
base_tree_level = 8
polydeg = 1

# Initial condition parameters for distribution function
epsilon = 0.01
k = 0.5

# SLURM job parameters
threads = 1
nnodes = 1
time = '16:00:00'
memory_request = '16G'

for M in Moments:
    command = f"julia examples/Vlasov_TwoStream.jl {M} {closure} {T_end} {base_tree_level} {polydeg} {epsilon} {k}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/VlasovPoisson/TwoStreamInstability/slurm_output_M{M}_closure{closure}_T_end{T_end}_epsilon{epsilon}_k{k}_p{polydeg}_level{base_tree_level}.out"
    )
