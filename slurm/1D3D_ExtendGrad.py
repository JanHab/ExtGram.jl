from clusterrun import createsbatch
import os

# Solver and simulation parameters
M = 4
closure = "Grad"
Kn = 1.0
source = "relaxation_source"
T_end = 0.3
base_tree_level = 10
polydeg = 1

# Initial condition parameters for the Riemann problem
rho_L_vec = [4.0, 7.0, 15.0]
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

# angles
angles = [ # maximizing angles
    3.14159, 1.5708,
    0.684719, 4.71239,
    2.03444, 1.5708,
    2.18628, 0.886077
]

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for rho_L in rho_L_vec:
    # 1D1D
    command = f"julia examples/Shock_1D1D.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L1} {theta_L} {rho_R} {v_R1} {theta_R} {x_lower} {x_upper}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/ShockTube/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L1}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R1}_theta_R{theta_R}.out"
    )


    # 1D3D
    command = f"julia examples/Shock_1D3D.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L1} {v_L2} {v_L3} {theta_L} {rho_R} {v_R1} {v_R2} {v_R3} {theta_R} {x_lower} {x_upper} {angles[0]} {angles[1]} {angles[2]} {angles[3]} {angles[4]} {angles[5]} {angles[6]} {angles[7]} {0}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/1D3D/Moments/slurm_output_1D3D_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L1}_{v_L2}_{v_L3}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R1}_{v_R2}_{v_R3}_theta_R{theta_R}.out"
    )
