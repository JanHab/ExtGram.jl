from clusterrun import createsbatch
import os

# Solver and simulation parameters
# M = 4
Moments = [4, 8, 12] #![4, 6, 8, 10, 12]
closure = "ExtGram"
Kn = 1.0
source = "relaxation_source"
T_end = 0.3
base_tree_level = 10
polydeg = 1
chi_vector = [-1.0, 0.0, 1.0, 1.5, "optimal", 2.0, 3.0]

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

for M in Moments:
    for chi_value in chi_vector:
        command = f"julia examples/GaugeInvestigation.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {chi_value} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/ShockTube/GaugeInvestigation/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_chi{chi_value}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
        )