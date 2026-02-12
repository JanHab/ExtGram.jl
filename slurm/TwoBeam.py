from clusterrun import createsbatch
import os

# Solver and simulation parameters
M_vector = [4, 5, 8, 9]
closure_vec = ["Gram", "ExtGram"]#!["Gram", "ExtGram", "Grad"]
Knudsen = [0.1, 1.0, 10.0] # to run
sources = ["relaxation_source", "relaxation_source", "zero_source"] # to run
T_end = 0.3
base_tree_level = 10
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_L = 1.0
v_L = 0.5
theta_L = 1.0
rho_R = 1.0
v_R = -0.5
theta_R = 1.0

x_lower = -2.0
x_upper = 2.0

# BGK
N = 500
c_l = -6.0
c_u = 6.0
x_left_BGK = -5.0
x_right_BGK = 5.0
base_tree_level_BGK = 10

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
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

    # BGK
    if source == "zero_source": # semi-analytic solution available, no need to run BGK
        continue
    command = f"julia examples/RiemannBGK.jl {N} {c_l} {c_u} {Kn} {source} {T_end} {base_tree_level_BGK} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_left_BGK} {x_right_BGK}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/Riemann1D/BGK/slurm_output_N{N}_c_l{c_l}_c_u{c_u}_Kn{Kn}_source{source}_T{T_end}_base_tree_level{base_tree_level_BGK}_polydeg{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}_x_left{x_left_BGK}_x_right{x_right_BGK}.out"
        )
