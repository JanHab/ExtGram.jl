from clusterrun import createsbatch
import os

# Solver and simulation parameters
Moments = [7, 9] #![4, 5, 6, 7, 8, 9]
closure = "ExtGram"
Knudsen = [0.2, 1.0]
source = "relaxation_source"

angle_name = "Max"

# SLURM job parameters
threads = 1
nnodes = 1
time = '02:00:00'
memory_request = '16G'

for i in range(len(Knudsen)):
    Kn = Knudsen[i]
    for M in Moments:
        time_ = time
        if closure == "Grad":
            time_ = '04:30:00'
        if closure == "ExtGram" and M == 7:
            time_ = '06:30:00'
        if closure == "ExtGram" and M == 9:
            time_ = '20:00:00'
        command = f"julia examples/Shock_1D3VSlab.jl {M} {closure} {Kn} {source} {angle_name}"
        # os.system(command)
        createsbatch(
            command, 
            nproc=threads, nnodes=nnodes, 
            time=time_, mem=memory_request, 
            output_file=f"out/1D3V/1D3V_angles/slurm_output_1D3V_M{M}_closure{closure}_Kn{Kn}_source{source}_anglepair{angle_name}.out"
        )

# closure = "Grad"
# Kn = 1.0
# time_grad = '04:30:00'
# for i in range(len(Knudsen)):
#     Kn = Knudsen[i]
#     for M in Moments:
#         command = f"julia examples/Shock_1D3VSlab.jl {M} {closure} {Kn} {source} {angle_name}"
#         # os.system(command)
#         createsbatch(
#             command, 
#             nproc=threads, nnodes=nnodes, 
#             time=time_grad, mem=memory_request, 
#             output_file=f"out/1D3V/1D3V_angles/slurm_output_1D3V_M{M}_closure{closure}_Kn{Kn}_source{source}_anglepair{angle_name}.out"
#             )

# One bad run for comparison 
# angle_name = "Bad"
# M = 4
# Kn = 1.0

# for closure in ["Grad"] #!["ExtGram", "Grad"]:
#     command = f"julia examples/Shock_1D3VSlab.jl {M} {closure} {Kn} {source} {angle_name}"
#     # os.system(command)
#     createsbatch(
#         command, 
#         nproc=threads, nnodes=nnodes, 
#         time=time, mem=memory_request, 
#         output_file=f"out/1D3V/1D3V_angles/slurm_output_1D3V_M{M}_closure{closure}_Kn{Kn}_source{source}_anglepair{angle_name}.out"
#         )
