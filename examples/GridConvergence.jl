using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = parse(Int, ARGS[1])    # number of moments
closure = ARGS[2] # flag for closure: "Gram", "ExtGram", "Grad"
Kn = parse(Float64, ARGS[3])  # Knudsen number
source_string = ARGS[4]
source = source_string == "relaxation_source" ? relaxation_source : zero_source # default to zero_source if not relaxation_source
T_end = parse(Float64, ARGS[5])

base_tree_level = parse(Int, ARGS[6]) # e.g. 8
polydeg = parse(Int, ARGS[7]) # e.g. 3

ρ_0 = parse(Float64, ARGS[8])
v_0 = parse(Float64, ARGS[9])
θ_0 = parse(Float64, ARGS[10])
ϵ = parse(Float64, ARGS[11])
k = parse(Float64, ARGS[12])

# Domain and discretization parameters
x_lower = parse(Float64, ARGS[13]); x_upper = parse(Float64, ARGS[14])
domain = (x_lower, x_upper)

equations = GramianMomentEquations1D(M, Kn, closure)

surface_flux = flux_lax_friedrichs
volume_flux = flux_central
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
basis = LobattoLegendreBasis(polydeg)
solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh(
    (domain[1],), (domain[2],), 
    initial_refinement_level=base_tree_level, 
    n_cells_max=100_000, 
    periodicity=true
)

initial_condition = InitialConditionsLandauDamping(
    ρ_0,
    ϵ,
    v_0,
    θ_0,
    k,
    equations
)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    source_terms=source
)

tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)


callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.99,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    time_interval = 500, # save at 500 time intervals
    name="GridConvergencePeriodic/moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name used for output
)

#= solve =#
sol = solve(
    ode, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    save_everystep=true,
    callback = callbacks,
);

# summary_callback()