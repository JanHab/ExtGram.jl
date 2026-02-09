from clusterrun import createsbatch
import os

# Solver and simulation parameters
M = 4
closure = "ExtGram"
Kn = 1.0
source = "relaxation_source"
T_end = 0.3
base_tree_level = 7#!10
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_L = 7.0
v_L1 = 0.0
v_L2 = 0.0
v_L3 = 0.0
theta_L = 1.0
rho_R = 1.0
v_R1 = 0.0
v_R2 = 0.0
v_R3 = 0.0
theta_R = 1.0

x_lower = -2.0
x_upper = 2.0

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

# command = f"julia examples/RiemannMoments1D3D.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L1} {v_L2} {v_L3} {theta_L} {rho_R} {v_R1} {v_R2} {v_R3} {theta_R} {x_lower} {x_upper}"
# os.system(command)
# createsbatch(
#     command, 
#     nproc=threads, nnodes=nnodes, 
#     time=time, mem=memory_request, 
#     output_file=f"out/Riemann1D/Moments/slurm_output_1D3D_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L1}_{v_L2}_{v_L3}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R1}_{v_R2}_{v_R3}_theta_R{theta_R}.out"
# )



# Checkin on angles
angles1 = [ # maximizing angles
    (3.14159, 1.5708),
    (0.684719, 4.71239),
    (2.03444, 1.5708),
    (2.18628, 0.886077)
]
angles2 = [ # maximizing angles / 2
    (3.14159/2, 1.5708/2),
    (0.684719/2, 4.71239/2),
    (2.03444/2, 1.5708/2),
    (2.18628/2, 0.886077/2)
]
angles3 = [
    (2.4, 0.3),
    (0.3, 3.314159),
    (1.3, 1.5),
    (2.9, 0.7)
]
angles4 = [
    (2.8, 1.1),
    (0.9, 4.3),
    (2.5, 1.1),
    (2.7, 0.6)
]
angles5 = [
    (2.3, 3.7),
    (1.3, 2.7),
    (2.9, 0.1),
    (0.7, 0.1)
]

angle_pairs = [angles1, angles2, angles3, angles4, angles5]

for i, angle_pair in enumerate(angle_pairs):
    command = f"julia examples/RiemannMoments1D3D_angles.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L1} {v_L2} {v_L3} {theta_L} {rho_R} {v_R1} {v_R2} {v_R3} {theta_R} {x_lower} {x_upper} {angle_pair} {i}"
    os.system(command)
    # createsbatch(
    #     command, 
    #     nproc=threads, nnodes=nnodes, 
    #     time=time, mem=memory_request, 
    #     output_file=f"out/Riemann1D/Moments/slurm_output_1D3D_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L1}_{v_L2}_{v_L3}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R1}_{v_R2}_{v_R3}_theta_R{theta_R}_anglepair{i}.out"
    # )