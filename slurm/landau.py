from clusterrun import createsbatch
import os

# Solver and simulation parameters
T_end = 25.0
base_tree_level = 8
polydeg = 1

# Initial condition parameters for distribution function
rho_0 = 1.0
v_0 = 0.0
theta_0 = 1.0
epsilon = 0.01
k = 0.5

# SLURM job parameters
threads = 1
nnodes = 1
time = '12:00:00'
memory_request = '16G'

# Varying closures for M=6
# Moments = [6, 7]
# closures = [
#     "Gram",
#     "ExtGram",
#     "Grad"
# ]
# for M in Moments:
#     for closure in closures:
#         command = f"julia examples/LandauDamping.jl {M} {closure} {T_end} {base_tree_level} {polydeg} {rho_0} {v_0} {theta_0} {epsilon} {k}"
#         # os.system(command)
#         createsbatch(
#             command, 
#             nproc=threads, nnodes=nnodes, 
#             time=time, mem=memory_request, 
#             output_file=f"out/VlasovPoisson/LandauDamping/slurm_output_M{M}_closure{closure}_T_end{T_end}_rho_0{rho_0}_v_0{v_0}_theta_0{theta_0}_epsilon{epsilon}_k{k}_p{polydeg}_level{base_tree_level}.out"
#         )



# Varying M with fixed closure "ExtGram"
Moments = [12, 13]#![5, 9, 13, 25]#![4, 5, 8, 9, 12, 24]
closure = 'ExtGram'

for M in Moments:
    command = f"julia examples/LandauDamping.jl {M} {closure} {T_end} {base_tree_level} {polydeg} {rho_0} {v_0} {theta_0} {epsilon} {k}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/VlasovPoisson/LandauDamping/slurm_output_M{M}_closure{closure}_T_end{T_end}_rho_0{rho_0}_v_0{v_0}_theta_0{theta_0}_epsilon{epsilon}_k{k}_p{polydeg}_level{base_tree_level}.out"
    )