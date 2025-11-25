from clusterrun import createsbatch

threads = 4
M_vector = [3] #! number of moments
closure = "Gram"
T_end = 25.0
base_tree_level = 8
polydeg = 3

for M in M_vector:
    command = f"julia examples/TwoStreamInstability_clusterrun.jl --threads={threads} --M={M} --closure={closure} --T_end={T_end} --base_tree_level={base_tree_level} --polydeg={polydeg}"
    createsbatch(command, nproc=threads, nnodes=1, time='00:15:00')