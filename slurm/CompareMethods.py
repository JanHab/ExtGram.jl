from clusterrun import createsbatch
import os

# Solver and simulation parameters
M = 8
# "ExtGram" already calculated from GridConvergence.py
# DVM already calculated from DVM.py
closure_vec = ["Gram", "ExtGram", "Grad"]
Knudsen = [0.1, 1.0, 10.0, 1.0] # last one doesn't matter, as zero_source
sources = ["relaxation_source", "relaxation_source", "relaxation_source", "zero_source"]
T_end = 0.3
base_tree_level = 8
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

# DVM
N = 1_000
c_l = -6.0
c_u = 6.0
x_left_DVM = -5.0
x_right_DVM = 5.0
base_tree_level_DVM = 8

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
    for closure in closure_vec:
        command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/ShockTube/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
        )
        if closure == "ExtGram":
            M_odd = 9
            command_odd = f"julia examples/RiemannMoments.jl {M_odd} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
            # os.system(command_odd)
            createsbatch(
                command_odd, 
                nproc=threads, nnodes=nnodes, 
                time=time, mem=memory_request, 
                output_file=f"out/ShockTube/Moments/slurm_output_M{M_odd}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
            )


    # DVM
    if source == "zero_source": # semi-analytic solution available, no need to run DVM
        continue
    command = f"julia examples/RiemannDVM.jl {N} {c_l} {c_u} {Kn} {source} {T_end} {base_tree_level_DVM} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_left_DVM} {x_right_DVM}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/ShockTube/DVM/slurm_output_N{N}_c_l{c_l}_c_u{c_u}_Kn{Kn}_source{source}_T{T_end}_base_tree_level{base_tree_level_DVM}_polydeg{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}_x_left{x_left_DVM}_x_right{x_right_DVM}.out"
        )



# Repeat with increased shock strength
Knudsen = [0.1, 1.0]
sources = ["relaxation_source", "relaxation_source"]

# Increase the shock
rho_L = 30.0


for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    source = sources[i]
    for closure in closure_vec:
        command = f"julia examples/RiemannMoments.jl {M} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time, mem=memory_request, 
            output_file=f"out/ShockTube/Moments/slurm_output_M{M}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
        )
        if closure == "ExtGram":
            M_odd = 9
            command_odd = f"julia examples/RiemannMoments.jl {M_odd} {closure} {Kn} {source} {T_end} {base_tree_level} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_lower} {x_upper}"
            # os.system(command_odd)
            createsbatch(
                command_odd, 
                nproc=threads, nnodes=nnodes, 
                time=time, mem=memory_request, 
                output_file=f"out/ShockTube/Moments/slurm_output_M{M_odd}_closure{closure}_Kn{Kn}_source{source}_T{T_end}_level{base_tree_level}_p{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}.out"
            )


    # DVM
    if source == "zero_source": # semi-analytic solution available, no need to run DVM
        continue
    command = f"julia examples/RiemannDVM.jl {N} {c_l} {c_u} {Kn} {source} {T_end} {base_tree_level_DVM} {polydeg} {rho_L} {v_L} {theta_L} {rho_R} {v_R} {theta_R} {x_left_DVM} {x_right_DVM}"
    # os.system(command)
    createsbatch(
        command, 
        nproc=threads, nnodes=nnodes, 
        time=time, mem=memory_request, 
        output_file=f"out/ShockTube/DVM/slurm_output_N{N}_c_l{c_l}_c_u{c_u}_Kn{Kn}_source{source}_T{T_end}_base_tree_level{base_tree_level_DVM}_polydeg{polydeg}_rho_L{rho_L}_v_L{v_L}_theta_L{theta_L}_rho_R{rho_R}_v_R{v_R}_theta_R{theta_R}_x_left{x_left_DVM}_x_right{x_right_DVM}.out"
        )