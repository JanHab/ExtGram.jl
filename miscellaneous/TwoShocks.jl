using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = 4    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.5
source = zero_source #!relaxation_source

# Domain and discretization parameters
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)
polydeg = 2  # polynomial degree
base_tree_level = 8  # initial mesh refinement level
surface_flux = flux_lax_friedrichs
volume_flux = flux_central

# Setting up everything
equations = GramianMomentEquations1D(M, Kn, closure)
initial_condition = InitialConditionsTwoShocks(
    Maxwellian(1.0, 0.0, 1.0), # Density, velocity, temperature
    Maxwellian(7.0, 0.0, 1.0), # Shock in density, but not velocity, temperature initially
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
# #= crashes for p > 1, i.e. this does not help at all
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    alpha_max = 1.0,#*0.5,#0.1, #  α_max = 1.0 seems natural -> corresponds to pure first order FV (Gassner paper)
    alpha_min = 0.01,#0.01,
    alpha_smooth = true, #* false, # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[3]#! *u[5]
) # ? Seems to restrict to M>=4
# `custom variable for smoothness detection?`

volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=false)

boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition), x_pos = BoundaryConditionDirichlet(initial_condition))

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
# savefig(p1, "out/Riemann1D/ρ_v_p.pdf")

# # Store primitive variables in CSV file
# CSV.write(
#     "out/Riemann1D/ρ_v_p.csv",
#     Tables.columntable((x=x, rho=ρ, v=v, p=p))
# )

# Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
n_plots = 5
p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
display(p2)
# savefig(p2, "out/Riemann1D/λ_max.pdf")

# Plot total variation in space over time
p3, TV_t = TVD_space(sol, M)
display(p3)
# savefig(p3, "out/Riemann1D/TV_space.pdf")
