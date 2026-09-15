
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LinearAlgebra, Statistics

closure = "ExtGram"
Kn = 0.5    # irrelevant
M = 4
theta_Max = [
    [3.14159, 0.684719, 2.03444, 2.18628], # M=4
    [3.14155, 2.57765, 2.28452, 0.684719, 1.95839, 2.03444], # M=6
    [3.14159, 0.490883, 0.729713, 1.91063, 0.955239, 0.857072, 2.57765, 1.1832, 2.03444] # M=8
]
phi_Max = [
    [1.5708, 4.71239, 1.5708, 0.886077], # M=4
    [1.5708, 1.5708, 1.5708, 5.27633, 1.5708, 2.13474], # M=6
    [1.5708, 4.71239, 4.71199, 1.5708, 4.71282, 4.22151, 1.07991, 4.22151, 0.841069] # M=8
]

theta_Arc = [
    [0., 0.589437, 1.57448, 2.02376], # M=4
    [0., 0.384521, 0.927303, 0.135409, 1.57322, 1.55268], # M=6
    [3.14159, 2.85885, 0.674741, 0.0587125, 2.03432, 2.27237, 1.56898, 1.57232, 1.78567] # M=8
]

phi_Arc = [
    [1.57079, 4.71239, 1.5708, 6.22436], # M=4
    [1.5708, 4.71239, 4.71239, 6.13099, 1.5708, 3.35395], # M=6
    [1.5708, 1.5708, 1.5708, 2.84783, 1.5708, 5.20013, 4.71239, 1.90366, 0.0886607] # M=8    
]

theta_Det = [
    [1.78445, 1.28569, 0.00017959, 0.701871], # M=4
    [2.41051, 1.83873, 2.68454, 0.836022, 0.069125, 0.411515], # M=6
    [1.96918, 0.934956, 1.90402, 0.284836, 0.708879, 0.18032, 0.363484, 1.7756, 2.4393] # M=8
]

phi_Det = [
    [0.790252, 1.50861, 1.36399, 2.03546], # M=4
    [1.68016, 2.1397, 5.38404, 3.71373, 5.01612, 4.59139], # M=6
    [5.45757, 3.62961, 4.55925, 1.24255, 4.07992, 5.22654, 2.21233, 4.65093, 0.740643] # M=8
]

# @testset "Compare Slab to Mathematica Reference" begin

ρ1 = 0.4
v1 = 0.0
θ1 = 0.6
ρ2 = 0.6
θ2 = 0.6
θ3 = 0.6

# E[Z^n] for Z ~ N(0,1): (n-1)!! for even n, 0 for odd n
std_normal_moment(n) = iseven(n) ? factorial(n) ÷ (2^(n ÷ 2) * factorial(n ÷ 2)) : 0

# Raw moment E[X^n] of X ~ N(μ, σ²), i.e. Mathematica's Moment[NormalDistribution[μ, Sqrt[σ2]], n]
moment1D(μ, σ2, n) = sum(binomial(n, k) * μ^(n - k) * σ2^(k ÷ 2) * std_normal_moment(k) for k in 0:2:n)

# Mathematica's Moment3D[{i, j, k}, v2]
moment3D((i, j, k), v2) =
    ρ1 * moment1D(v1, θ1, i) * moment1D(0.0, θ3, j) * moment1D(0.0, θ3, k) +
    ρ2 * moment1D(v2, θ2, i) * moment1D(0.0, θ3, j) * moment1D(0.0, θ3, k)

# Full slab state vector in the ordering used by GramianMomentEquations1D3V
moments3D(v2, M, slab_geometry) = [moment3D(Tuple(ix), v2) for ix in U_t_index(M, slab_geometry)]






# Mathematica's Sqrt[Mean[((approx - exact)/exact)^2]]: RMS of the element-wise relative error
relative_rms_error(approx, exact) = sqrt(mean(((approx .- exact) ./ exact) .^ 2))





# v2_vec = range(0.1, 4.0, length=20)
v2_vec = [0.1, 0.3, 0.8, 1.5, 2.5, 3.0]
angles = ["Max", "Det", "Arc"]
# !angles = ["Max", "Det"]
for M in [4, 6, 8]
    relative_RMS_error = Float64[]
    relL2 = Float64[]

    for angle in angles
        if angle == "Max"
            theta = theta_Max[M ÷ 2 - 1];
            phi = phi_Max[M ÷ 2 - 1];
        elseif angle == "Arc"
            theta = theta_Arc[M ÷ 2 - 1];
            phi = phi_Arc[M ÷ 2 - 1];
        elseif angle == "Det"
            theta = theta_Det[M ÷ 2 - 1];
            phi = phi_Det[M ÷ 2 - 1];
        end
        equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=true, theta=theta, phi=phi);

        for v2 in v2_vec
            u = moments3D(v2, M, true);
            closure_eval = ExtGram.closure_moments(u, equations);
            closure_moments_ground_truth = [moment3D(Tuple(ix), v2) for ix in equations._closure_index];
            error_rms = relative_rms_error(closure_eval, closure_moments_ground_truth);
            error_l2 = norm(closure_eval - closure_moments_ground_truth, 2) / norm(closure_moments_ground_truth, 2);
            push!(relative_RMS_error, error_rms);
            push!(relL2, error_l2);
        end
    end

    plt = plot(xlabel="v2", ylabel="Relative RMS Error", title="Closure Error vs v2 for M=$M", yscale=:log10)
    plot!(plt, v2_vec, relative_RMS_error[1:length(v2_vec)], label="Max")
    plot!(plt, v2_vec, relative_RMS_error[length(v2_vec)+1:2*length(v2_vec)], label="Det")
    plot!(plt, v2_vec, relative_RMS_error[2*length(v2_vec)+1:end], label="Arc")
    # display(plt)
    savefig(plt, "closure_error_M_$M.png")
end


# Testing my moment assembly vs the mathematica implementation
test = [1., 0.6, 0, 0, 1.2, 0., 0., 0.6, 0, 0.6, 1.68, 0, 0, 0.36, 0., 0.36, 0, 0, 0, 0, 3.84, 0., 0., 0.72, 0, 0.72, 0., 0., 0., 0., 1.08, 0, 0.36, 0, 1.08, 7.44, 0, 0, 1.008, 0., 1.008, 0, 0, 0, 0, 0.648, 0., 0.216, 0., 0.648, 0, 0, 0, 0, 0, 0, 18.96, 0., 0., 2.304, 0, 2.304, 0., 0., 0., 0., 1.296, 0, 0.432, 0, 1.296, 0., 0., 0., 0., 0., 0., 3.24, 0, 0.648, 0, 0.648, 0, 3.24]
testing_my_implelentaiton = moments3D(1.0, 6, false)
println(isapprox(test, testing_my_implelentaiton, rtol=1e-4))


test2 = [1., 1.8, 0, 0, 6., 0., 0., 0.6, 0, 0.6, 19.44, 0, 0, 1.08, 0., 1.08, 0, 0, 0, 0, 69.12, 0., 0., 3.6, 0, 3.6, 0., 0., 0., 0., 1.08, 0, 0.36, 0, 1.08, 252.72, 0, 0, 11.664, 0., 11.664, 0, 0, 0, 0, 1.944, 0., 0.648, 0., 1.944, 0, 0, 0, 0, 0, 0, 965.52, 0., 0., 41.472, 0, 41.472, 0., 0., 0., 0., 6.48, 0, 2.16, 0, 6.48, 0., 0., 0., 0., 0., 0., 3.24, 0, 0.648, 0, 0.648, 0, 3.24, 3802.46, 0, 0, 151.632, 0., 151.632, 0, 0, 0, 0, 20.9952, 0., 6.9984, 0., 20.9952, 0, 0, 0, 0, 0, 0, 5.832, 0., 1.1664, 0., 1.1664, 0., 5.832, 0, 0, 0, 0, 0, 0, 0, 0, 15462.6, 0., 0., 579.312, 0, 579.312, 0., 0., 0., 0., 74.6496, 0, 24.8832, 0, 74.6496, 0., 0., 0., 0., 0., 0., 19.44, 0, 3.888, 0, 3.888, 0, 19.44, 0., 0., 0., 0., 0., 0., 0., 0., 13.608, 0, 1.944, 0, 1.1664, 0, 1.944, 0, 13.608]
testing_my_implelentaiton2 = moments3D(3.0, 8, false)
println(isapprox(test2, testing_my_implelentaiton2, rtol=1e-4))