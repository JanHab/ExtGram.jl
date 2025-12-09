using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = 4 #parse(Int, ARGS[1])
closure = "ExtGram" #ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = 1.0 #parse(Float64, ARGS[3])
# source_string = ARGS[4]
source = zero_source #relaxation_source
T_end = 0.3 #parse(Float64, ARGS[5])
base_tree_level = 2 #parse(Int, ARGS[6]) # e.g. 8
polydeg = 1 #parse(Int, ARGS[7]) # e.g. 3
ρ_L = 7.0 #parse(Float64, ARGS[8]) # 7.0
v_L1 = 0.0 #parse(Float64, ARGS[9]) # 0.0
v_L2 = 0.0
v_L3 = 0.0
θ_L = 1.0 #parse(Float64, ARGS[10]) # 1.0
ρ_R = 1.0 #parse(Float64, ARGS[11]) # 1.0
v_R1 = 0.0 #parse(Float64, ARGS[12]) # 0.0
v_R2 = 0.0
v_R3 = 0.0
θ_R = 1.0 #parse(Float64, ARGS[13]) # 1.0


# Fixed settings
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)



# Setting up everything
equations = GramianMomentEquations1D3D(M, Kn, closure)
# Mp1 = 35
# f_left = Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L)
# left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
initial_condition = InitialConditionsShockTube1D3D(
    Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L), # Density, velocity, temperature
    Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R), # Shock in density, but not velocity, temperature initially
    M,
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
# #= crashes for p > 1, i.e. this does not help at all
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    alpha_max = 0.5, #!1.0,#*0.5,#0.1, #  α_max = 1.0 seems natural -> corresponds to pure first order FV (Gassner paper)
    alpha_min = 0.001, #!0.01,#0.01,
    alpha_smooth = true, #* false, # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[5]#! *u[5]
) # ? Seems to restrict to M>=4

surface_flux = flux_lax_friedrichs
volume_flux = flux_central
volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh(
    (domain[1],), (domain[2],), 
    initial_refinement_level=base_tree_level, 
    n_cells_max=10_000, 
    periodicity=false
)

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

# callbacks, summary_callback = callbacksGramianMomentEquations(
#     semi, tspan, basis; 
#     cfl = 0.9,          # Maximum cfl number
#     plot_interval = 20,  # plot every 20 steps
#     name="gram_solution",
# )

cfl = 0.99
time_interval = 20
name = "gram_solution_1D3D"

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=cfl)

save_solution = SaveTriangulationCallback(
    time_interval=tspan[2]/time_interval,
    save_initial_solution=true,
    file_format="tsv",
    append_solution=true,
    solution_variables = cons2cons,
    clear_out_dir=false,
    name=name,
    info="basis = $(Base.typename(typeof(basis)).wrapper)"
)

callbacks = CallbackSet(
    alive_callback,
    stepsize_callback,
    # plot_callback,
    save_solution,
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
# x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)
# display(p1)
# savefig(p1, "out/Riemann1D/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_rho_v_p.pdf")

# # Store primitive variables in CSV file
# CSV.write(
#     "out/Riemann1D/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_rho_v_p.csv",
#     Tables.columntable((x=x, rho=ρ, v=v, p=p))
# )

# # Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
# n_plots = 5
# p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
# display(p2)
# savefig(p2, "out/Riemann1D/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_lambda_max.pdf")