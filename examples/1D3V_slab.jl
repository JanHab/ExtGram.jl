
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables


closure = "ExtGram"
Kn = 0.5    # irrelevant
M = 4
theta_Max = [
    # [3.14159, 0.684719, 2.03444, 2.18628], # M=4
    [3.14159, 0.616578, 1.10715, 2.18628, ], # M=4, mirrored
    # [3.14155, 2.57765, 2.28452, 0.684719, 1.95839, 2.03444], # M=6
    [3.14155, 2.80135, 0.857072, 0.684719, 1.1832, 2.0944], # M=6, mirrored
    # [3.14159, 0.490883, 0.729713, 1.91063, 0.955239, 0.857072, 2.57765, 1.1832, 2.03444] # M=8
    [3.14159, 0.390413, 0.0000975558, 2.57765, 3.14159, 0.785398, 1.91063, 1.10711, 1.10715] # M=8, mirrored
]
phi_Max = [
    # [1.5708, 4.71239, 1.5708, 0.886077], # M=4
    [1.5708, 4.39298, 4.71239, 0.886077], # M=4, mirrored
    # [1.5708, 1.5708, 1.5708, 5.27633, 1.5708, 2.13474], # M=6
    [1.5708, 2.02967, 4.71239, 5.27633, 4.71239, 2.28452], # M=6, mirrored
    # [1.5708, 4.71239, 4.71199, 1.5708, 4.71282, 4.22151, 1.07991, 4.22151, 0.841069] # M=8
    [1.5708, 5.01795, 3.98292, 1.07991, 2.52616, 5.32787, 1.5708, 5.44221, 3.98266]
]

theta_Arc = [
    # [0., 0.589437, 1.57448, 2.02376], # M=4
    [5.85555e-10, 3.14159, 3.14159, 1.57601], # M=4, mirrored
    # [0., 0.384521, 0.927303, 0.135409, 1.57322, 1.55268], # M=6
    [0., 3.14159, 3.14159, 0.593484, 3.14159, 1.57352], # M=6, mirrored
    # [3.14159, 2.85885, 0.674741, 0.0587125, 2.03432, 2.27237, 1.56898, 1.57232, 1.78567] # M=8
    [0., 0., 0., 0.386383, 3.14159, 0.785703, 0., 1.56657, 1.56918], # M=8, mirrored
]

phi_Arc = [
    # [1.57079, 4.71239, 1.5708, 6.22436], # M=4
    [4.71239, 5.30183, 0.00368397, 2.35619], # M=4, mirrored
    # [1.5708, 4.71239, 4.71239, 6.13099, 1.5708, 3.35395], # M=6
    [1.5708, 1.18628, 3.78509, 1.06088, 0.0024264, 3.92699], # M=6, mirrored
    # [1.5708, 1.5708, 1.5708, 2.84783, 1.5708, 5.20013, 4.71239, 1.90366, 0.0886607] # M=8    
    [1.5708, 4.99513, 5.38713, 4.35201, 3.60512, 2.18642, 0.00181229, 5.94573, 2.35619] # M=8, mirrored
]

theta_Det = [
    # [1.78445, 1.28569, 0.00017959, 0.701871], # M=4
    [1.08562, 0.969517, 0.0841159, 0.140558], # M=4, mirrored
    # [2.41051, 1.83873, 2.68454, 0.836022, 0.069125, 0.411515], # M=6
    [0.791178, 2.71441, 1.19527, 0.7375, 2.69368, 2.97235], # M=6, mirrored
    # [1.96918, 0.934956, 1.90402, 0.284836, 0.708879, 0.18032, 0.363484, 1.7756, 2.4393] # M=8
    [0.787033, 2.64384, 0.0612564, 2.99817, 0.313342, 0.516051, 2.44407, 0.0226567, 0.99256], # M=8, mirrored
]

phi_Det = [
    # [0.790252, 1.50861, 1.36399, 2.03546], # M=4
    [1.57834, 2.11307, 4.62546, 4.0341], # M=4, mirrored
    # [1.68016, 2.1397, 5.38404, 3.71373, 5.01612, 4.59139], # M=6
    [1.49439, 5.46196, 0.779535, 0.457368, 4.84181, 1.53241], # M=6, mirrored
    # [5.45757, 3.62961, 4.55925, 1.24255, 4.07992, 5.22654, 2.21233, 4.65093, 0.740643] # M=8
    [5.32313, 4.40688, 5.47647, 5.13619, 0.599894, 5.59516, 2.9642, 1.72031, 2.16574], # M=8, mirrored
]

# @testset "Compare Slab to Mathematica Reference" begin

ρ1 = 0.4
v1 = 0.0
θ1 = 0.6
ρ2 = 0.6
θ2 = 0.6
θ3 = 0.6

# function moments3D(v2, M)
#     f1 = Maxwellian1D3V(ρ1, (v1, 0.0, 0.0), θ1)
#     f2 = Maxwellian1D3V(ρ2, (v2, 0.0, 0.0), θ2)
    
#     full = ExtGram.multi_index_list(M)
#     equations = GramianMomentEquations1D3V(M, 1.0, "ExtGram", slab_geometry=true)#, theta=theta_Max[1], phi=phi_Max[1])# Only relevant is the slab_geometry=true to getht the correct positons
#     position_in_full = Dict(ix => i for (i, ix) in enumerate(full))
#     gather = [position_in_full[ix] for ix in equations._U_t_index]
#     u1 = convective_moments_1D3V(M, f1.ρ, f1.v, f1.θ)[gather]
#     u2 = convective_moments_1D3V(M, f2.ρ, f2.v, f2.θ)[gather]
#     return u1 .+ u2
# end

# E[Z^n] for Z ~ N(0,1): (n-1)!! for even n, 0 for odd n
std_normal_moment(n) = iseven(n) ? factorial(n) ÷ (2^(n ÷ 2) * factorial(n ÷ 2)) : 0

# Raw moment E[X^n] of X ~ N(μ, σ²), i.e. Mathematica's Moment[NormalDistribution[μ, Sqrt[σ2]], n]
moment1D(μ, σ2, n) = sum(binomial(n, k) * μ^(n - k) * σ2^(k ÷ 2) * std_normal_moment(k) for k in 0:2:n)

# Mathematica's Moment3D[{i, j, k}, v2]
moment3D((i, j, k), v2) =
    ρ1 * moment1D(v1, θ1, i) * moment1D(0.0, θ3, j) * moment1D(0.0, θ3, k) +
    ρ2 * moment1D(v2, θ2, i) * moment1D(0.0, θ3, j) * moment1D(0.0, θ3, k)

# Full slab state vector in the ordering used by GramianMomentEquations1D3V
moments3D(v2, M) = [moment3D(Tuple(ix), v2) for ix in U_t_index(M, true)]


# # v2=1.0 and M=4
v2 = 6.0
M = 4
u = moments3D(v2, M)
# u = [
#     1.0, # M=0
#     0.6, # M=1
#     1.2, 0.6, # M=2
#     1.68, 0.36, # M=3
#     3.84, 0.72, 1.08, 0.36 # M=4
# ]

closure_moments_mathematica = [
    [7.48604, 1.0158, 0.648279, 0.211783],
    [7.48604, 1.0158, 0.648, 0.213438],
    [7.48604, 1.00752, 0.649537, 0.215253]
]


for angle_name in ["Max", "Arc", "Det"]
    if angle_name == "Max"
        theta = theta_Max[M ÷ 2 - 1]
        phi = phi_Max[M ÷ 2 - 1]
    elseif angle_name == "Arc"
        theta = theta_Arc[M ÷ 2 - 1]
        phi = phi_Arc[M ÷ 2 - 1]
    elseif angle_name == "Det"
        theta = theta_Det[M ÷ 2 - 1]
        phi = phi_Det[M ÷ 2 - 1]
    else
        error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", or \"Det\".")
    end

    equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=true, theta=theta, phi=phi)
    closure_eval = ExtGram.closure_moments(u, equations)


    println(isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2))
    closure_moments_ground_truth = [moment3D(Tuple(ix), v2) for ix in equations._closure_index]
    println(isapprox(closure_eval, closure_moments_ground_truth, rtol=1e-2))

    println("Moments predicted: $closure_eval")
    println("Moments ground truth: $closure_moments_ground_truth")

    # closure_moments_ground_truth = moments3D(v2, M+1)[ExtGram.index_1d(M+1)]
    # print(isapprox(closure_eval, closure_moments_ground_truth, rtol=1e-2))
end

v2 = 6.0
M = 6
u = moments3D(v2, M)
# u = [
#     1.0, # M=0
#     0.6, # M=1
#     1.2, 0.6, # M=2
#     1.68, 0.36, # M=3
#     3.84, 0.72, 1.08, 0.36, # M=4
#     7.44, 1.008, 0.648, 0.216, # M=5
#     18.96, 2.304, 1.296, 0.432, 3.24, 0.648, # M=6
# ]

closure_moments_mathematica = [
    [44.5385, 4.4743, 1.81966, 0.612977, 1.94435, 0.39299],
    [44.5385, 4.47436, 1.81979, 0.619747, 1.944, 0.388748],
    [44.5385, 0.445648, 1.39957, 1.69309, 2.20915, 0.278901]
]

# for angle_name in ["Max", "Arc", "Det"]
angle_name = "Max"
    if angle_name == "Max"
        theta = theta_Max[M ÷ 2 - 1]
        phi = phi_Max[M ÷ 2 - 1]
    elseif angle_name == "Arc"
        theta = theta_Arc[M ÷ 2 - 1]
        phi = phi_Arc[M ÷ 2 - 1]
    elseif angle_name == "Det"
        theta = theta_Det[M ÷ 2 - 1]
        phi = phi_Det[M ÷ 2 - 1]
    else
        error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", or \"Det\".")
    end

    equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=true, theta=theta, phi=phi)
    closure_eval = ExtGram.closure_moments(u, equations)

    println(isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2))
    closure_moments_ground_truth = [moment3D(Tuple(ix), v2) for ix in equations._closure_index]
    println(isapprox(closure_eval, closure_moments_ground_truth, rtol=1e-2))

    println("Moments predicted: $closure_eval")
    println("Moments ground truth: $closure_moments_ground_truth")

    # closure_moments_ground_truth = moments3D(v2, M+1)[ExtGram.index_1d(M+1)]
    # print(isapprox(closure_eval, closure_moments_ground_truth, rtol=1e-2))
end


# v2=1.0 and M=8
M = 8
u = moments3D(v2, M)
# u = [
#     1.0, # M=0
#     0.6, # M=1
#     1.2, 0.6, # M=2
#     1.68, 0.36, # M=3
#     3.84, 0.72, 1.08, 0.36, # M=4
#     7.44, 1.008, 0.648, 0.216, # M=5
#     18.96, 2.304, 1.296, 0.432, 3.24, 0.648, # M=6
#     44.448, 4.464, 1.8144, 0.6048, 1.944, 0.3888, # M=7
#     124.08, 11.376, 4.1472, 1.3824, 3.888, 0.7776, 13.608, 1.944, 1.1664 # M=8
# ]

closure_moments_mathematica = [
    [332.214, 26.691, 8.04271, 2.68905, 5.44984, 1.09466, 8.16533, 1.24152, 0.558735],
    [332.214, 26.6909, 8.04265, 2.68115, 5.45014, 1.07082, 8.1648, 1.16638, 0.699903],
    [332.976, 25.7487, 10.5176, 2.69411, 5.68976, 0.864099, 8.09875, 0.356056, 0.849342]
]

for angle_name in ["Max", "Arc", "Det"]
    if angle_name == "Max"
        theta = theta_Max[M ÷ 2 - 1]
        phi = phi_Max[M ÷ 2 - 1]
    elseif angle_name == "Arc"
        theta = theta_Arc[M ÷ 2 - 1]
        phi = phi_Arc[M ÷ 2 - 1]
    elseif angle_name == "Det"
        theta = theta_Det[M ÷ 2 - 1]
        phi = phi_Det[M ÷ 2 - 1]
    else
        error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", or \"Det\".")
    end

    equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=true, theta=theta, phi=phi)
    closure_eval = ExtGram.closure_moments(u, equations)

    println(isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2))
end