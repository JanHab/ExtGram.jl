from clusterrun import createsbatch
import os

# Solver and simulation parameters
M = 4
closure = "ExtGram"
T_end = 25.0
base_tree_level = 6
polydeg = 1 #!3

# Initial condition parameters for distribution function
rho_0 = 1.0
v_0 = 0.0
theta_0 = 1.0
epsilon = 0.01
k = 0.5

# SLURM job parameters
threads = 1
nnodes = 1
time = '04:00:00'
memory_request = '16G'

command = f"julia examples/LandauDamping.jl {M} {closure} {T_end} {base_tree_level} {polydeg} {rho_0} {v_0} {theta_0} {epsilon} {k}"
os.system(command)
# createsbatch(
#     command, 
#     nproc=threads, nnodes=nnodes, 
#     time=time, mem=memory_request, 
#     output_file=f"out/VlasovPoisson/LandauDamping/slurm_output_M{M}_closure{closure}_T_end{T_end}_rho_0{rho_0}_v_0{v_0}_theta_0{theta_0}_epsilon{epsilon}_k{k}_p{polydeg}_level{base_tree_level}.out"
# )