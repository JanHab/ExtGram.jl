using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = 4
closure = "ExtGram" # String # "Gram", "ExtGram" or "Grad"
Kn = 1.0
source = zero_source
ρ_L = 7.0
v_L1 = 1.0
v_L2 = 0.0
v_L3 = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R1 = 1.0
v_R2 = 0.0
v_R3 = 0.0
θ_R = 1.0

# 1D1D IC
f_left = Maxwellian(ρ_L, v_L1, θ_L)
left_1D1D = convective_moments(f_left, Val(M+1))
f_right = Maxwellian(ρ_R, v_R1, θ_R)
right_1D1D = convective_moments(f_right, Val(M+1))

# 1D3D IC
# slab_indices = Int[]
# index_start = 0
# for i in 0:M
#     append!(slab_indices, HyQMOM.index_1d(i) .+ index_start) # Offset by shell start
#     index_start += size(HyQMOM.index(i), 1)
# end
slab_indices = get_valid_indices(M)
f_left = Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L)
left_1D3D = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
left_1D3D = left_1D3D[slab_indices]
f_right = Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R)
right_1D3D = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)
right_1D3D = right_1D3D[slab_indices]

println("Initial Moments (Left Side, 1D1D / 1D3D):")
println("First Order Moments: ", left_1D1D[1], " / ", left_1D3D[1])
println("Second Order Moments: ", left_1D1D[2], " / ", left_1D3D[2])
println("Third Order Moments: ", left_1D1D[3], " / ", left_1D3D[3:4])
println("Fourth Order Moments: ", left_1D1D[4], " / ", left_1D3D[5:6])
println("Fifth Order Moments: ", left_1D1D[5], " / ", left_1D3D[7:10])

println("\nInitial Moments (Right Side, 1D1D / 1D3D):")
println("First Order Moments: ", right_1D1D[1], " / ", right_1D3D[1])
println("Second Order Moments: ", right_1D1D[2], " / ", right_1D3D[2])
println("Third Order Moments: ", right_1D1D[3], " / ", right_1D3D[3:4])
println("Fourth Order Moments: ", right_1D1D[4], " / ", right_1D3D[5:6])
println("Fifth Order Moments: ", right_1D1D[5], " / ", right_1D3D[7:10])

# Closure test
equations_1D1D = GramianMomentEquations1D(M, Kn, closure)
val_1D1D_left = HyQMOM.closure(left_1D1D, equations_1D1D)
val_1D1D_right = HyQMOM.closure(right_1D1D, equations_1D1D)

equations_1D3D = GramianMomentEquations1D3D(M, Kn, closure)
val_1D3D_left = HyQMOM.closure_transform(left_1D3D, equations_1D3D)
val_1D3D_right = HyQMOM.closure_transform(right_1D3D, equations_1D3D)

println("\nClosure Results (Left Side, 1D1D / 1D3D):")
println("Closure Moment: ", val_1D1D_left, " / ", val_1D3D_left)
println("\nClosure Results (Right Side, 1D1D / 1D3D):")
println("Closure Moment: ", val_1D1D_right, " / ", val_1D3D_right)
