from clusterrun import createsbatch
import os
import numpy as np

# Solver and simulation parameters
M_vector = [4, 5]
closure = "ExtGram"
Kn = 1.0
source = "relaxation_source"
T_end = 0.3
base_tree_level_vec = [6, 7, 8, 9, 10, 11, 12, 14]
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_0 = 1.0
v_0 = 0.0
theta_0 = 1.0
epsilon = 0.01
k = 0.5

x_lower = 0.0
x_upper = 4.0 * np.pi

# SLURM job parameters
threads = 1
nnodes = 1
time = '16:00:00'
memory_request = '32G'

for M in M_vector:
    for base_tree_level in base_tree_level_vec:
        command = f"julia examples/GridConvergence.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_0} {v_0} {theta_0} {epsilon} {k} {x_lower} {x_upper}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/GridConvergencePeriodic/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_0{rho_0}_v_0{v_0}_theta_0{theta_0}_epsilon{epsilon}_k{k}.out"
        )