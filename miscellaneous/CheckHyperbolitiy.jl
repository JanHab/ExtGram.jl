using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

M = 4
closure = "ExtGram" # "Gram", "ExtGram" or "Grad"
Kn = 1.0 # not important here, just needed for seting up the equations

equations = GramianMomentEquations1D3D(M, Kn, closure)

# IC
v_R2_vec = LinRange(0.0, 3.0, 10)
# ρ_L_vec = LinRange(0.4, 0.6, 3)
# v_L1_vec = LinRange(0.0, 1.0, 10)
# v_L3_vec = LinRange(0.0, 1.0, 10)
for v2 in v_R2_vec
    ρ_L = 0.4
    v_L1 = 0.0
    v_L2 = 0.1 #!0.0
    v_L3 = 0.0
    θ_L = 0.6

    ρ_R = 0.6
    v_R1 = v2#!1.5
    v_R2 = 0.0
    v_R3 = 0.0
    θ_R = 0.6

    f_left = Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L)
    moments_init = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
    f_right = Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R)
    moments_init_right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)
    moments_init = moments_init .+ moments_init_right
    # moments_init = [1., 0.9, 0, 0, 1.95, 0., 0., 0.6, 0, 0.6, 3.645, 0, 0, 0.54, 0., 0.54, 0, 0, 0, 0, 8.9775, 0., 0., 1.17, 0, 1.17, 0., 0., 0., 0., 1.08, 0, 0.36, 0, 1.08]

    # Reduce to slab geometry
    slab_indices = []
    index_start = 0
    for i in 0:M
        append!(slab_indices, index_1d(i) .+ index_start) # Offset by shell start
        index_start += size(index(i), 1)
    end
    # moments_init = [i in slab_indices ? moments_init[i] : 0.0 for i in 1:length(moments_init)]
    moments_init = moments_init[slab_indices]

    ExtGram.closure_transform(moments_init, equations)

    ####################################
    # Check Hyperbolicity of the system #
    ####################################
    flux_jac = ExtGram.flux_jacobian(moments_init, equations)
    eigenvals = eigvals(flux_jac)
    # println("Eigenvalues of the flux Jacobian:")
    # println(eigenvals)
    # println("Are all eigenvalues real? ", all(isreal, eigenvals))
    real_parts = real.(eigenvals)
    imag_parts = imag.(eigenvals)
    # scatter(real_parts, imag_parts, xlabel="Real Part", ylabel="Imaginary Part", title="Eigenvalues of Flux Jacobian", legend=false)
    maximum_imag = maximum(abs.(imag_parts))
    if maximum_imag < 1e-10
        println("The system is hyperbolic (all eigenvalues are real).")
    else
        println("The system is not hyperbolic (some eigenvalues are complex).")
        println("Maximum imaginary part of eigenvalues: ", maximum_imag)
        println("All Eigenvalues:")
        println(eigenvals)
    end
end

display(flux_jac)

real_parts = real.(eigenvals)
imag_parts = imag.(eigenvals)
scatter(real_parts, imag_parts, xlabel="Real Part", ylabel="Imaginary Part", title="Eigenvalues of Flux Jacobian", legend=false)
maximum_imag = maximum(abs.(imag_parts))