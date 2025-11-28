from clusterrun import createsbatch
import os

threads = 1 #!16
# M_vector = [14] #![5, 7, 9, 11, 13]#![3, 5, 7, 9, 11, 13] #! number of moments
# closures = ["ExtGram"]#!["Gram", "ExtGram"]
M_vector = [4, 8, 14, 5, 9, 15]
closures = ["Gram", "ExtGram"]
T_end = 25.0
base_tree_level = 8
polydeg = 1 #!3

nnodes = 1
time = '04:00:00'
memory_request = '16G'

for closure in closures:
    for M in M_vector:
        # command = f"julia examples/TwoStreamInstability_clusterrun.jl {M} {closure} {T_end} {base_tree_level} {polydeg} --threads={threads}"
        command = f"julia examples/TwoStreamInstability_clusterrun.jl {M} {closure} {T_end} {base_tree_level} {polydeg}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/VlasovPoisson/TwoStreamInstability/slurm_output_closure{closure}_T{T_end}_M{M}_p{polydeg}_level{base_tree_level}.out"
        )