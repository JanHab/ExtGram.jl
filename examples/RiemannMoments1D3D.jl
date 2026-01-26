using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = 4 #parse(Int, ARGS[1])
closure = "ExtGram" #ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = 1.0 #parse(Float64, ARGS[3])
# source_string = ARGS[4]
source_string = "zero_source"#!"relaxation_source"
source = zero_source #!relaxation_source #zero_source #relaxation_source
T_end = 25 #parse(Float64, ARGS[5])
base_tree_level = 8 #!2 #parse(Int, ARGS[6]) # e.g. 8
polydeg = 1 #parse(Int, ARGS[7]) # e.g. 3

# Riemann
# ρ_L = 7.0 #parse(Float64, ARGS[8]) # 7.0
# v_L1 = 0.0 #!1.0 #!1.5 #parse(Float64, ARGS[9]) # 0.0
# v_L2 = 0.0
# v_L3 = 0.0
# θ_L = 1.0 #parse(Float64, ARGS[10]) # 1.0
# ρ_R = 1.0 #parse(Float64, ARGS[11]) # 1.0
# v_R1 = 0.0 #!1.0 #parse(Float64, ARGS[12]) # 0.0
# v_R2 = 0.0
# v_R3 = 0.0
# θ_R = 1.0 #parse(Float64, ARGS[13]) # 1.0

# x_lower = -2.0; x_upper = 2.0

# Rankine-Hugoniot
Ma = 1.4#!3.8

ρ_L = 1.0
ρ_R = (4*Ma^2) / (Ma^2 + 3)

v_L1 = sqrt(5/3) * Ma
v_L2 = 0.0
v_L3 = 0.0
v_R1 = sqrt(5/3) * (Ma^2 + 3.0) / (4.0 * Ma)
v_R2 = 0.0
v_R3 = 0.0

θ_L = 1.0
θ_R = ( (5*Ma^2 - 1) * (Ma^2 + 3) ) / (16 * Ma^2)

# x_lower = -30; x_upper = 30
x_lower = -50; x_upper = 50




domain = (x_lower, x_upper)



# Setting up everything
equations = GramianMomentEquations1D3D(M, Kn, closure)
# Mp1 = 35
# f_left = Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L)
# left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
# f_right = Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R)
# right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)
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
    variable = (u, eqns)->u[1]*u[3]#! *u[5]
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
# name = "gram_solution_1D3D"
name = "1D3D/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L1$(v_L1)_v_L2$(v_L2)_v_L3$(v_L3)_v_R1$(v_R1)_v_R2$(v_R2)_v_R3$(v_R3)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)"

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=cfl)

plot_interval = 20  # plot every 20 steps
plot_callback = VisualizationCallback(
    semi;
    interval=plot_interval,
    solution_variables=cons2cons,
    plot_data_creator=PlotData1D,
    plot_creator=Trixi.show_plot,
)

save_solution_cons = SaveTriangulationCallback(
    time_interval=tspan[2]/time_interval,
    save_initial_solution=true,
    file_format="tsv",
    append_solution=true,
    solution_variables = cons2cons,
    clear_out_dir=false,
    name=name * "_cons",
    info="basis = $(Base.typename(typeof(basis)).wrapper)"
)
save_solution_prim = SaveTriangulationCallback(
    time_interval=tspan[2]/time_interval,
    save_initial_solution=true,
    file_format="tsv",
    append_solution=true,
    solution_variables = cons2prim,
    clear_out_dir=false,
    name=name * "_prim",
    info="basis = $(Base.typename(typeof(basis)).wrapper)"
)

callbacks = CallbackSet(
    alive_callback,
    stepsize_callback,
    plot_callback,
    #! save_solution_cons,
    #! save_solution_prim,
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
    # saveat = range(tspan[1], tspan[2], length=100)
);

# # summary_callback()

# Post Processing
n_equations = 10
u_final = sol.u[end] #!ode.u0 #!sol.u[end]
L = length(u_final)
Nloc = L ÷ (n_equations)
Fmat = reshape(u_final, n_equations, Nloc)

u = []
for i in 1:n_equations
    push!(u, Fmat[i, :])
end

x = range(x_lower, x_upper, length=Nloc)
ρ = vec(u[1, :])[1]
v = vec(u[2, :])[1] ./ ρ
θ = (vec(u[3, :])[1] + vec(u[4, :])[1] + vec(u[4, :])[1]) ./ (3.0 .* ρ)
# p = ρ .* θ
# θ = p ./ ρ
# Plot primitive variables at final time
plot(x, (ρ .- ρ_L) ./ (ρ_R - ρ_L), label="ρ", lw=2, xlim=(-10, 10))
plot(x, (v .- v_R1) ./ (v_L1 - v_R1), label="v", lw=2, xlim=(-10, 10))
# plot(x, (p .- θ_L .* ρ_L) ./ (θ_R .* ρ_R - θ_L .* ρ_L), label="p", lw=2, xlim=(-10, 10))
plot(x, (θ .- θ_L) ./ (θ_R - θ_L), label="θ", lw=2, xlim=(-10, 10))

pl = plot(xlim=(-10, 10), title="T=$(T_end)", size=(500,500), yticks=0:0.1:1);
plot!(pl, x, (ρ .- ρ_L) ./ (ρ_R - ρ_L), label="ρ", lw=2);
plot!(pl, x, (v .- v_R1) ./ (v_L1 - v_R1), label="v", lw=2);
# plot!(pl, x, (p .- θ_L .* ρ_L) ./ (θ_R .* ρ_R - θ_L .* ρ_L), label="p", lw=2)
plot!(pl, x, (θ .- θ_L) ./ (θ_R - θ_L), label="θ", lw=2);
display(pl)

# Plot first moments
pl = plot();
plot!(
    pl,
    xlabel="x",
    ylabel="U_i",
);
x = range(x_lower, x_upper, length=Nloc)

plot!(
    pl,
    x, u[1, :], 
    label="U000",
    color=:blue,
);

plot!(
    pl,
    x, u[2, :], 
    label="U100",
    color=:red,
);

plot!(
    pl,
    x, u[3, :], 
    label="U200",
    color=:green,
);

# plot!(
#     pl,
#     x, u[4, :], 
#     label="U020",
#     color=:orange,
# )
plot!(
    pl,
    x, u[5, :], 
    label="U300",
    color=:purple,
);

plot!(
    pl,
    x, u[7, :], 
    label="U400",
    color=:brown,
);

display(pl)
