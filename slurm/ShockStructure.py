from clusterrun import createsbatch
import os
import numpy as np

# Solver and simulation parameters
Moments = [4]#!, 8]
closure_vec = ["Gram"]#!, "ExtGram", "Grad"]
Kn = 1.0
source = "relaxation_source"
T_end = 50.0 #!25.0
base_tree_level = 8 #!10
polydeg = 1

# Initial condition parameters for the RH problem
Mach_numbers = [1.4]#!, 2.0]

x_left = -20.0
x_right = 100.0

# BGK
N = 500
c_l = -10.0
c_u = 10.0

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for Ma in Mach_numbers:
    rho_L = 1.0
    rho_R = (2*Ma**2) / (Ma**2 + 1)
    v_L = np.sqrt(3) * Ma
    v_R = np.sqrt(3) / 2 * (Ma**2 + 1) / Ma
    theta_L = 1.0
    theta_R = (3*Ma**2 - 1) * (Ma**2 + 1) / (4 * Ma**2)
    # Moment Methods
    for closure in closure_vec:
        for M in Moments:
            command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_left} {x_right}"
            os.system(command)
    #         createsbatch(
    #             command, 
    #             nproc=threads, nnodes=nnodes, 
    #             time=time, mem=memory_request, 
    #             output_file=f"out/Riemann1D/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
    #         )

    # # BGK
    # command = f"julia examples/RiemannBGK.jl {N} {c_l} {c_u} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_left} {x_right}"
    # # os.system(command)
    # createsbatch(
    #     command, 
    #     nproc=threads, nnodes=nnodes, 
    #     time=time, mem=memory_request, 
    #     output_file=f"out/Riemann1D/BGK/slurm_output_N{N}_c_l{c_l}_c_u{c_u}_Kn{Kn}_source{source}_T{T_end}_base_tree_level{base_tree_level}_polydeg{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}_x_left{x_left}_x_right{x_right}.out"
    #     )
