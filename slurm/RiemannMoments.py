from clusterrun import createsbatch
import os

# Solver and simulation parameters
M = 24
closure = "ExtGram"
Kn = 1.0
source = "zero_source"#!"relaxation_source"
T_end = 0.3
base_tree_level = 6
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_L = 7.0
v_L = 0.0
theta_L = 1.0
rho_R = 1.0
v_R = 0.0
theta_R = 1.0

# SLURM job parameters
threads = 1
nnodes = 1
time = '04:00:00'
memory_request = '16G'

command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R}"
os.system(command)
# createsbatch(
#     command, 
#     nproc=threads, nnodes=nnodes, 
#     time=time, mem=memory_request, 
#     output_file=f"out/Riemann1D/Moments/relaxation/slurm_output__M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
# )