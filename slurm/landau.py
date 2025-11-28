from clusterrun import createsbatch
import os

threads = 1 #!16
M_vector = [25]#![5, 15]
closures = ["Gram", "ExtGram"]#!, "Grad"]
T_end = 25.0
base_tree_level = 6
polydeg = 3

nnodes = 1
time = '04:00:00'
memory_request = '16G'

for closure in closures:
    for M in M_vector:
        command = f"julia examples/LandauDamping_clusterrun.jl {M} {closure} {T_end} {base_tree_level} {polydeg}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/VlasovPoisson/LandauDamping/slurm_output_closure{closure}_T{T_end}_M{M}_p{polydeg}_level{base_tree_level}.out"
        )