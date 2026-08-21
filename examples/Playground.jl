
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays




M = 4
U_t_index_ = U_t_index(M, true)
U_x_index_ = U_x_index(M, true)


V = [@SVector[A[i,1], A[i,2], A[i,3], A[i,4]] for i in axes(A,1)]

# * The angles -> t.b.d.
# theta_Nmax = [3.14159, 0.463648, 0.684719, 3.14159, 2.25552, 0.61548, 2.03444, 2.18628, 0.0000534309, 1.5708, 1.5708, 1.5708]
# phi_Nmax = [1.5708, 4.71239, 4.71239, 0.886077, 4.71239, 4.02767, 1.5708, 0.886077, 3.60524, 4.71084, 4.02767, 3.60524]

# theta_Arc = [0.0, 0.244979, 0.588003, 3.14159, 0.982794, 2.67795, 1.81578, 1.10715, 0.0, 1.5708, 1.5708, 1.5708]
# phi_Arc = [1.5708, 4.71239, 1.5708, 0.981359, 4.71239, 2.30207, 4.71239, 2.30207, 6.2795, 4.71239, 5.30183, 6.2795]

# theta_det = [0.973951, 0.399601, 0.587204, 2.84498, 2.08866, 1.77566, 0.33008, 2.32745, 1.78321, 2.78056, 0.206615, 2.83077]
# phi_det = [5.00082, 1.67262, 1.40439, 0.223503, 5.77741, 1.60421, 2.47911, 2.27102, 1.08665, 4.49793, 3.89513, 2.46093]

# Slab-geometry angles
theta_Nmax = [3.14159, 0.684719, 2.034440, 2.186280]
phi_Nmax = [1.5708, 4.71239, 1.570800, 0.886077]

theta, phi = theta_Nmax, phi_Nmax 

M = 4
Kn = 1.0
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
domain = (-2.0, 2.0)
base_tree_level = 4 # ! Increase later
source = zero_source
T_end = 0.025 #! Increase later 0.3
cfl = 0.99 # ? Decrease later? 0.45
time_interval = 10 # ? Increase later? 10
name = "out"

equations = GramianMomentEquations1D3V(M, Kn, "ExtGram", theta=theta, phi=phi, slab_geometry=true)

f_left = Maxwellian1D3D(7.0, (0.0, 0.0, 0.0), 1.0)
f_right = Maxwellian1D3D(1.0, (0.0, 0.0, 0.0), 1.0)

initial_condition = InitialConditionsShockTube1D3V(
    f_left, # Density, velocity, temperature
    f_right, # Shock in density, but not velocity, temperature initially
    M,
    equations
)

basis = LobattoLegendreBasis(polydeg)

# shock capturing
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    variable = (u, eqns)->u[1]*u[5] # * Placeholder for the moment
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
stepsize_callback = StepsizeCallback(cfl=cfl)

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
    alive_callback, summary_callback, stepsize_callback,
    save_solution_cons, save_solution_prim
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


# sol_end_cons = sol.u[end]
# sol_end_prim = cons2prim(sol_end_cons, equations)
# x = vec(semi.cache.elements.node_coordinates)

# ρ = sol_end_cons[1]
# v = sol_end_cons[2] / ρ
# P_200 = sol_end_prim[5]
# P_020 = sol_end_prim[8]
# P_002 = sol_end_prim[10]
# θ = 1 / (3 * ρ) * (P_200 + P_020 + P_002 - ρ * v^2)



# p = plot(sol.t, [ρ, v, θ], label=["Density" "Velocity" "Temperature"], xlabel="Time", ylabel="Values", title="Shock Tube Evolution", lw=2)