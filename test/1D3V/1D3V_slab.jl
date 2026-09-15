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

# v2=1.0 and M=4
u = [
    1.0, # M=0
    0.6, # M=1
    1.2, 0.6, # M=2
    1.68, 0.36, # M=3
    3.84, 0.72, 1.08, 0.36 # M=4
]

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

    @test isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2)
end

# v2=1.0 and M=6
M = 6
u = [
    1.0, # M=0
    0.6, # M=1
    1.2, 0.6, # M=2
    1.68, 0.36, # M=3
    3.84, 0.72, 1.08, 0.36, # M=4
    7.44, 1.008, 0.648, 0.216, # M=5
    18.96, 2.304, 1.296, 0.432, 3.24, 0.648, # M=6
]

closure_moments_mathematica = [
    [44.5385, 4.4743, 1.81966, 0.612977, 1.94435, 0.39299],
    [44.5385, 4.47436, 1.81979, 0.619747, 1.944, 0.388748],
    [44.5385, 0.445648, 1.39957, 1.69309, 2.20915, 0.278901]
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

    @test isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2)
end


# v2=1.0 and M=8
M = 8
u = [
    1.0, # M=0
    0.6, # M=1
    1.2, 0.6, # M=2
    1.68, 0.36, # M=3
    3.84, 0.72, 1.08, 0.36, # M=4
    7.44, 1.008, 0.648, 0.216, # M=5
    18.96, 2.304, 1.296, 0.432, 3.24, 0.648, # M=6
    44.448, 4.464, 1.8144, 0.6048, 1.944, 0.3888, # M=7
    124.08, 11.376, 4.1472, 1.3824, 3.888, 0.7776, 13.608, 1.944, 1.1664 # M=8
]

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

    @test isapprox(closure_eval, closure_moments_mathematica[M ÷ 2 - 1], rtol=1e-2)
end