using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = 4    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.3
source = relaxation_source
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)

# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DRiemann(
    M, Kn, closure,
    Maxwellian(7.0, 0.0, 1.0), # Density, velocity, temperature
    Maxwellian(1.0, 0.0, 1.0);
    domain = domain,
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
    name="Riemann1D/gram_solution",
)



#= solve =#
sol = solve(
    ode, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    callback = callbacks,
    saveat = range(tspan[1], tspan[2], length=100)
);

summary_callback()


# Post Processing
x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)
display(p1)
savefig(p1, "out/Riemann1D/ρ_v_p.pdf")

# Store primitive variables in CSV file
CSV.write(
    "out/Riemann1D/ρ_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)

# Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
n_plots = 5
p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
display(p2)
savefig(p2, "out/Riemann1D/λ_max.pdf")

# Plot total variation in space over time
p3, TV_t = TVD_space(sol, M)
display(p3)
savefig(p3, "out/Riemann1D/TV_space.pdf")
