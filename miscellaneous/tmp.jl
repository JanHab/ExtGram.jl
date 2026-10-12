
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi , OrdinaryDiffEq, Plots #, CSV, Tables, LinearAlgebra

using LinearAlgebra, StaticArrays

# Access arguments by index
M = 4
closure = "ExtGram"
Kn = 1.0
source_string = "relaxation_source"
source = source_string == "relaxation_source" ? relaxation_source : zero_source # default to zero_source if not relaxation_source
T_end = 0.3
base_tree_level = 10
polydeg = 1
ρ_L = 7.0
v_L = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R = 0.0
θ_R = 1.0
x_lower = -2.0
x_upper = 2.0
volume_flux = flux_central
surface_flux = flux_lax_friedrichs

# Fixed settings
domain = (x_lower, x_upper)

equations = GramianMomentEquations1D(M, Kn, closure)
initial_condition = InitialConditionsShockTube(
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R), # Density, velocity, temperature
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    variable = (u, eqns)->u[1]*u[3]
)

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

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=0.45)

plot_callback = VisualizationCallback(
    semi;
    interval=20,
    solution_variables=cons2cons,
    plot_data_creator=PlotData1D,
    plot_creator=Trixi.show_plot
)

callbacks = CallbackSet(
    alive_callback, summary_callback, stepsize_callback, plot_callback
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