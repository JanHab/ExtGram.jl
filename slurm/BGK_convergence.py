from clusterrun import createsbatch
import os

threads = 1
N_vector = [10, 50, 100, 200]#!, 500, 1_000]
T_end = 0.2
base_tree_level = 8
polydeg = 1
Kn = 0.1

nnodes = 1
time = '04:00:00'
memory_request = '16G'

for N in N_vector:
    command = f"julia examples/Riemann1DBGK_clusterrun.jl {N} {T_end} {base_tree_level} {polydeg} {Kn}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/Riemann1D/slurm_output_Kn{Kn}_T{T_end}_N{N}_p{polydeg}_level{base_tree_level}.out"
        )
