using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
N = 50
c_l = -6.0
c_u = 6.0
Kn = 1.0 # Knudsen number
source = relaxation_source

domain = (-5.0, 5.0)
T_end = 0.5

basis, mesh, equations, initial_condition, solver, boundary_conditions = setupBGK1DRiemann(
    N, Kn,
    c_l, c_u, 
    Maxwellian(7.0, 0.0, 1.0), # Density, velocity, temperature
    Maxwellian(1.0, 0.0, 1.0);
    base_tree_level = 6,
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
    name="Riemann1D/bgk_solution",
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
# Plotting density, velocity, 
x = LinRange(domain[1], domain[2], length(sol.u[end]) ÷ N)
ρ, v, θ, p = ρ_v_θ_p_BGK(sol.u[end], equations)

p1 = plot_ρ_v_p_bgk(ρ, v, p, x; xlims=(-2.0, 2.0))
display(p1)
savefig(p1, "out/Riemann1D/bgk_ρ_v_p.pdf")

# Store primitive variables in CSV file
CSV.write(
    "out/Riemann1D/bgk_ρ_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)
