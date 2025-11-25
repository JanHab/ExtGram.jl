from clusterrun import createsbatch

threads = 16
M_vector = [4] #![3, 5, 7, 9, 11, 13] #! number of moments
closure = "Gram"
T_end = 25.0
base_tree_level = 8
polydeg = 3

nnodes = 1
time = '00:15:00'

for M in M_vector:
    command = f"julia examples/TwoStreamInstability_clusterrun.jl {M} {closure} {T_end} {base_tree_level} {polydeg} --threads={threads}"
    createsbatch(command, nproc=threads, nnodes=nnodes, time=time, outfile=f"out/VlasovPoisson/TwoStreamInstability/slurm_output_closure{closure}_T{T_end}_M{M}_p{polydeg}_level{base_tree_level}.out")