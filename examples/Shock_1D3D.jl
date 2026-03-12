"""
    Example file for a 1D-3D shock (tube or structure) test case with the moment method.

    ! Note: Only implemented for M=4 for now, but can be extended to higher M in the future.
"""

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = parse(Int, ARGS[1])
@assert M ==4 # !only M=4 is currently supported for 1D3D
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = parse(Float64, ARGS[3])
source_string = ARGS[4]
source = source_string == "relaxation_source" ? relaxation_source : zero_source 
T_end = parse(Float64, ARGS[5])
base_tree_level = parse(Int, ARGS[6]) # e.g. 10
polydeg = parse(Int, ARGS[7]) # e.g. 1

# Riemann
ρ_L = parse(Float64, ARGS[8]) # 7.0
v_L1 = parse(Float64, ARGS[9]) # 0.0
v_L2 = parse(Float64, ARGS[10]) # 0.0
v_L3 = parse(Float64, ARGS[11]) # 0.0
θ_L = parse(Float64, ARGS[12]) # 1.0
ρ_R = parse(Float64, ARGS[13]) # 1.0
v_R1 = parse(Float64, ARGS[14]) # 0.0
v_R2 = parse(Float64, ARGS[15]) # 0.0
v_R3 = parse(Float64, ARGS[16]) # 0.0
θ_R = parse(Float64, ARGS[17]) # 1.0


# Fixed settings
x_lower = parse(Float64, ARGS[18]) # -2.0
x_upper = parse(Float64, ARGS[19]) # 2.0
domain = (x_lower, x_upper)

angles = [ # maximizing angles
    (parse(Float64, ARGS[20]), parse(Float64, ARGS[21])),
    (parse(Float64, ARGS[22]), parse(Float64, ARGS[23])),
    (parse(Float64, ARGS[24]), parse(Float64, ARGS[25])),
    (parse(Float64, ARGS[26]), parse(Float64, ARGS[27]))
]

anglepairnumber = parse(Int, ARGS[28]) 



# Setting up everything
equations = GramianMomentEquations1D3D(M, Kn, closure, angles=angles)

initial_condition = InitialConditionsShockTube1D3D(
    Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L), # Density, velocity, temperature
    Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R), # Shock in density, but not velocity, temperature initially
    M,
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    alpha_max = 0.15, 
    alpha_min = 0.001, 
    alpha_smooth = true,  # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[5]
)

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

cfl = 0.99
time_interval = 20
name = "1D3D/1D3D_angles/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L1$(v_L1)_v_L2$(v_L2)_v_L3$(v_L3)_v_R1$(v_R1)_v_R2$(v_R2)_v_R3$(v_R3)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_anglepairnumber$(anglepairnumber)"

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
    save_solution_cons,
    save_solution_prim,
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
);

summary_callback()


# Eigenvalues
if polydeg == 1 # only implemented for polydeg=1 (linear basis functions)
    n_equations = 10
    u_final = sol.u[end]
    L = length(u_final)
    Nloc = L ÷ (n_equations)
    Fmat = reshape(u_final, n_equations, Nloc)
    x = range(x_lower, x_upper, length=Nloc)

    u = []
    for i in 1:n_equations
        push!(u, Fmat[i, :])
    end

    x_vector = [];
    λ_vector = [];
    u_vector = [];
    for x_index in 1:Nloc
        u = SVector{n_equations}(Fmat[:, x_index]...)
        jacobian = flux_jacobian(u, equations)
        λ = real.(eigen(jacobian).values)
        
        push!(x_vector, x[x_index])
        push!(λ_vector, λ)
        push!(u_vector, u)
    end

    CSV.write(
        "out/1D3D/1D3D_angles/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L1$(v_L1)_v_L2$(v_L2)_v_L3$(v_L3)_v_R1$(v_R1)_v_R2$(v_R2)_v_R3$(v_R3)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_anglepairnumber$(anglepairnumber)_eigenvalues.csv",
        Tables.columntable((x_vector=x_vector, eigenvalue=λ_vector, moments=u_vector))
    )
end