from clusterrun import createsbatch
import os

# Solver and simulation parameters
M_vector = [4, 5]#!, 6, 7, 8, 9, 10, 12]
closure_vec = ["Gram", "ExtGram", "Grad"]
Knudsen = [0.1, 1.0]#!, 10.0]
sources = ["relaxation_source", "relaxation_source"]#!, "zero_source"]
T_end = 0.3
base_tree_level = 10
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_L = 7.0
v_L = 0.0
theta_L = 1.0
rho_R = 1.0
v_R = 0.0
theta_R = 1.0

x_lower = -2.0
x_upper = 2.0

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
    for closure in closure_vec:
        for M in M_vector:
            command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
            # os.system(command)
            createsbatch(
                command, 
                nproc=threads, nnodes=nnodes, 
                time=time, mem=memory_request, 
                output_file=f"out/Riemann1D/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
            )



# # Increase M for Kn = 1.0
# M_vector = [6, 7, 8, 9, 10, 12]
# closure_vec = ["Gram", "ExtGram"]#!, "Grad"]
# Kn = 1.0
# source = "relaxation_source"

# for closure in closure_vec:
#     for M in M_vector:
#         command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
#         # os.system(command)
#         createsbatch(
#             command, 
#             nproc=threads, nnodes=nnodes, 
#             time=time, mem=memory_request, 
#             output_file=f"out/Riemann1D/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
#         )