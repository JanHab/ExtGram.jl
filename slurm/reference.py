from clusterrun import createsbatch
import os

# SLURM job parameters
threads = 1
nnodes = 1
time = '08:00:00'
memory_request = '16G'

# rho_L_vec = [7.0, 30.0]
# Kn_vec = [0.2, 1.0]

# nx_vec = [128, 256, 512]#!, 1024]
# nu_vec = [100, 200, 400]
# nx x nu: 512 x 400 (~40 min), 512 x 800 (~75 min), 1024 x 800 (~5 h)
nx = 512 
nu = 400
rho_L_vec = [7.0, 30.0]
Kn = 1.0

for rho_L in rho_L_vec:
    # 1D-3V DVM reference solution with KitBase.jl
    command = f"julia -t 4 examples/DVM_Reference.jl {Kn} {rho_L} {nx} {nu}"
    os.system(command)
