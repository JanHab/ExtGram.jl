using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = 4    # number of moments
extended = true # flag for extended gramian closure or "standard" closure
Kn = 1.0  # Knudsen number
T_end = 1.0#!2.5
source = vlasov_poisson_source
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)

# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DRiemann(
    M, Kn, extended,
    Maxwellian(1.0, 0.0, 1.0), # Density, velocity, temperature
    Maxwellian(1.0, 0.0, 1.0);
    domain = domain,
    base_tree_level=8
)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    boundary_conditions=boundary_conditions, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
)
# Add Vlasov-Poisson callback
callbacks = CallbackSet(callbacks, vlasov_poisson_callback(;M, domain))

#= solve =#
sol = solve(
    ode, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    callback = callbacks,
);

summary_callback()

# Access the energy history
energy_history = HyQMOM.ELECTRIC_FIELD.Energy

# You can then plot it or analyze it
# todo: fix to cfl (not fixed time-step)
time = LinRange(0, T_end, length(energy_history))
plot(
    time, 
    energy_history ./ energy_history[1], 
    xlabel="Time", ylabel="Normalized Electric Field Energy", 
    label="Energy",
    # yaxis=:log
)
