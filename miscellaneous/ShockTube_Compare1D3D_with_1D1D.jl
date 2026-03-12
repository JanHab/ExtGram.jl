"""
    Some routine to run the 1D-1D and 1D-3D Riemann problem with the same initial conditions and compare the results.
        - The same closure is used for both.
"""

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

##################################################################
###################### ## ShockTube Problem ######################
##################################################################
M = 4 
closure = "ExtGram" 
Kn = 1.0 
source = relaxation_source 
T_end = 0.3
base_tree_level = 8
polydeg = 1 
ρ_L = 7.0 
v_L1 = 0.0 
v_L2 = 0.0
v_L3 = 0.0
θ_L = 1.0 
ρ_R = 1.0 
v_R1 = 0.0 
v_R2 = 0.0
v_R3 = 0.0
θ_R = 1.0 


# Fixed settings
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)



# Setting up everything
equations_1D3D = GramianMomentEquations1D3D(M, Kn, closure)
equations_1D1D = GramianMomentEquations1D(M, Kn, closure)
initial_condition_1D3D = InitialConditionsShockTube1D3D(
    Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L), # Density, velocity, temperature
    Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R), # Shock in density, but not velocity, temperature initially
    M,
    equations_1D3D
)
initial_condition_1D1D = InitialConditionsShockTube(
    Maxwellian(ρ_L, v_L1, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R1, θ_R), # Shock in density, but not velocity, temperature initially
    equations_1D1D
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
# #= crashes for p > 1, i.e. this does not help at all
indicator_sc_1D3D = IndicatorHennemannGassner(
    equations_1D3D, basis,
    alpha_max = 0.5, # Rather conservative guess
    alpha_min = 0.001, 
    alpha_smooth = true, # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[5]
) 
indicator_sc_1D1D = IndicatorHennemannGassner(
    equations_1D1D, basis,
    alpha_max = 0.5, # Rather conservative guess
    alpha_min = 0.001, 
    alpha_smooth = true, # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[3]
) 

#= surface flux and volume flux =#
surface_flux = flux_lax_friedrichs
volume_flux = flux_central
volume_integral_1D3D = VolumeIntegralShockCapturingHG(
    indicator_sc_1D3D;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)
volume_integral_1D1D = VolumeIntegralShockCapturingHG(
    indicator_sc_1D1D;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

# = solver =#
solver_1D3D = DGSEM(basis, surface_flux, volume_integral_1D3D)
solver_1D1D = DGSEM(basis, surface_flux, volume_integral_1D1D)

# = mesh =#
mesh = TreeMesh(
    (domain[1],), (domain[2],), 
    initial_refinement_level=base_tree_level, 
    n_cells_max=10_000, 
    periodicity=false
)

# = boundary conditions =#
boundary_conditions_1D3D = (x_neg = BoundaryConditionDirichlet(initial_condition_1D3D), x_pos = BoundaryConditionDirichlet(initial_condition_1D3D))
boundary_conditions_1D1D = (x_neg = BoundaryConditionDirichlet(initial_condition_1D1D), x_pos = BoundaryConditionDirichlet(initial_condition_1D1D))

# = semidiscretization =#
semi_1D3D = SemidiscretizationHyperbolic(
    mesh, equations_1D3D, 
    initial_condition_1D3D, solver_1D3D, 
    boundary_conditions=boundary_conditions_1D3D, 
    source_terms=source
)
semi_1D1D = SemidiscretizationHyperbolic(
    mesh, equations_1D1D, 
    initial_condition_1D1D, solver_1D1D, 
    boundary_conditions=boundary_conditions_1D1D, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode_1D3D = semidiscretize(semi_1D3D, tspan)
ode_1D1D = semidiscretize(semi_1D1D, tspan)

cfl = 0.99
time_interval = 20

#= callbacks =#
alive_callback_1D3D = AliveCallback(analysis_interval=100)
summary_callback_1D3D = SummaryCallback()
stepsize_callback_1D3D = StepsizeCallback(cfl=cfl)

save_solution_1D3D = SaveTriangulationCallback(
    time_interval=tspan[2]/time_interval,
    save_initial_solution=true,
    file_format="tsv",
    append_solution=true,
    solution_variables = cons2cons,
    clear_out_dir=false,
    name="gram_solution_1D3D",
    info="basis = $(Base.typename(typeof(basis)).wrapper)"
)

callbacks_1D3D = CallbackSet(
    alive_callback_1D3D,
    stepsize_callback_1D3D,
    # plot_callback,
    save_solution_1D3D,
)
callbacks_1D1D, summary_callback_1D1D = callbacksGramianMomentEquations(
    semi_1D1D, tspan, basis; 
    cfl = 0.99,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
)

#= solve =#
sol_1D3D = solve(
    ode_1D3D, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()...,
    callback = callbacks_1D3D,
    # saveat = range(tspan[1], tspan[2], length=100)
);
sol_1D1D = solve(
    ode_1D1D, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    callback = callbacks_1D1D,
    # saveat = range(tspan[1], tspan[2], length=100)
);


# Post Processing
n_equations_1D3D = 10
u_final_1D3D = sol_1D3D.u[end]
u_final_1D1D = sol_1D1D.u[end]
L_1D3D = length(u_final_1D3D)
L_1D1D = length(u_final_1D1D)
Nloc = L_1D3D ÷ (n_equations_1D3D)
Fmat_1D3D = reshape(u_final_1D3D, n_equations_1D3D, Nloc)
Fmat_1D1D = reshape(u_final_1D1D, M+1, Nloc)

u_1D3D = []
for i in 1:n_equations_1D3D
    push!(u_1D3D, Fmat_1D3D[i, :])
end
u_1D1D = []
for i in 1:(M+1)
    push!(u_1D1D, Fmat_1D1D[i, :])
end

# Plot moments
x = range(x_lower, x_upper, length=Nloc)
pl = Vector{Any}(undef, M+1)
for i in 1:M+1
    pl[i] = scatter()
    plot!(
        pl[i],
        xlabel="x",
        ylabel="Moment $(i-1)",
    )
end

plot!(
    pl[1],
    x, u_1D3D[1], 
    label="1D3D: U000",
    color=:blue,
)
plot!(
    pl[1],
    x, u_1D1D[1], 
    label="1D1D: u0",
    color=:cyan,
)

pl[2] = plot()
plot!(
    pl[2],
    x, u_1D3D[2], 
    label="1D3D: U100",
    color=:red,
)
plot!(
    pl[2],
    x, u_1D1D[2], 
    label="1D1D: u1",
    color=:orange,
)

pl[3] = plot()
plot!(
    pl[3],
    x, u_1D3D[3], 
    label="1D3D: U200",
    color=:green,
)
plot!(
    pl[3],
    x, u_1D1D[3], 
    label="1D1D: u2",
    color=:lime,
)

pl[4] = plot()
plot!(
    pl[4],
    x, u_1D3D[5], 
    label="1D3D: U300",
    color=:purple,
)
plot!(
    pl[4],
    x, u_1D1D[4], 
    label="1D1D: u3",
    color=:magenta,
)

pl[5] = plot()
plot!(
    pl[5],
    x, u_1D3D[7], 
    label="1D3D: U400",
    color=:brown,
)
plot!(
    pl[5],
    x, u_1D1D[5], 
    label="1D1D: u4",
    color=:pink,
)

l = @layout  [grid(2,3)]# a{0.2w}]
plot(pl..., layout=l, size=(1_200, 600))
