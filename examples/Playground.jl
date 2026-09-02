
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays




# * The angles -> t.b.d.
# theta_Nmax = [3.14159, 0.463648, 0.684719, 3.14159, 2.25552, 0.61548, 2.03444, 2.18628, 0.0000534309, 1.5708, 1.5708, 1.5708]
# phi_Nmax = [1.5708, 4.71239, 4.71239, 0.886077, 4.71239, 4.02767, 1.5708, 0.886077, 3.60524, 4.71084, 4.02767, 3.60524]

# theta_Arc = [0.0, 0.244979, 0.588003, 3.14159, 0.982794, 2.67795, 1.81578, 1.10715, 0.0, 1.5708, 1.5708, 1.5708]
# phi_Arc = [1.5708, 4.71239, 1.5708, 0.981359, 4.71239, 2.30207, 4.71239, 2.30207, 6.2795, 4.71239, 5.30183, 6.2795]

# theta_det = [0.973951, 0.399601, 0.587204, 2.84498, 2.08866, 1.77566, 0.33008, 2.32745, 1.78321, 2.78056, 0.206615, 2.83077]
# phi_det = [5.00082, 1.67262, 1.40439, 0.223503, 5.77741, 1.60421, 2.47911, 2.27102, 1.08665, 4.49793, 3.89513, 2.46093]

# Slab-geometry angles
# theta_Nmax = [3.14159, 0.684719, 2.034440, 2.186280]
# phi_Nmax = [1.5708, 4.71239, 1.570800, 0.886077]

# theta, phi = theta_Nmax, phi_Nmax 

M = 4
Kn = 1.0
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
domain = (-2.0, 2.0)
base_tree_level = 10 # ! Increase later
source = zero_source # ? relaxation_source #zero_source
T_end = 0.3 #! Increase later 0.3
cfl = 0.99 # ? Decrease later? 0.45
time_interval = 10 # ? Increase later? 10
name = "out"

# ? theta = [0.042978, 2.86589, 1.76672, 1.62279, 0.599645, 2.28761, 2.84739, 0.25247, 0.194722, 2.28605, 0.00700916, 2.83612, 2.00972, 2.7489, 2.16911]
# ? phi = [3.77752, 3.82357, 1.1951, 5.16362, 1.54305, 5.36209, 4.51818, 2.59034, 2.02184, 3.98618, 4.83726, 5.44366, 3.69455, 3.42223, 4.85428]

theta21 = [2.84178, 2.25739, 2.83628, 2.27054, 1.51937, 1.55188, 2.50057, 0.396858, 1.92276, 1.27198, 1.01401, 0.739861, 1.09696, 0.332375, 2.21655, 1.58578, 1.62733, 0.779441, 0.126096, 3.09066, 2.97756]
phi21 = [4.38201, 5.52717, 1.89082, 3.76511, 4.05576, 1.45790, 2.99889, 4.45293, 3.41820, 4.91036, 2.76854, 4.70358, 4.48148, 3.88086, 4.72682, 2.42491, 1.96357, 1.95470, 1.57480, 5.33683, 1.10492]

equations = GramianMomentEquations1D3V(M, Kn, "ExtGram", theta=theta21, phi=phi21)
# equations = GramianMomentEquations1D3V(M, Kn, "ExtGram") #!, theta=theta, phi=phi, slab_geometry=true)

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


#= ---------- final state: ρ, v, p, θ over x ---------- =#
u_end = Array(Trixi.wrap_array(sol.u[end], semi))          # (nvars, nnodes, nelements)
x     = vec(Array(semi.cache.elements.node_coordinates))   # same (node, element) ordering

nvars = nvariables(equations)
prim  = vec([cons2prim(SVector{nvars}(u_end[:, i, e]), equations)
             for i in axes(u_end, 2), e in axes(u_end, 3)])

i200, i020, i002 = equations._pressure_index      # (5, 8, 10) for M = 4

ρ   = [p[1] for p in prim]
v_x   = [p[2] for p in prim]
v_y = [p[3] for p in prim]
v_z = [p[4] for p in prim]
θ   = [(p[i200] + p[i020] + p[i002]) / (3 * p[1]) for p in prim]
prs = ρ .* θ                                       # scalar pressure p = ρθ

perm = sortperm(x); xs = x[perm]
plt = plot(
    plot(xs, ρ[perm], ylabel="ρ", title="density"),
    plot(xs, [v_x[perm] v_y[perm] v_z[perm]], label=["v_x" "v_y" "v_z"], title="velocity", legend=:topright),
    plot(xs, prs[perm], ylabel="p", title="pressure"),
    plot(xs, θ[perm], ylabel="θ", title="temperature"),
    layout=(2,2), size=(900,600), xlabel="x", lw=2, legend=false,
    plot_title="1D-3V  M=$M  Kn=$Kn  t=$(round(sol.t[end], digits=3))"
)
display(plt)





M_vec = [4, 5, 6, 7, 8]
for M in M_vec
    name = "1D3V_full/M=$M"

    equations = GramianMomentEquations1D3V(M, Kn, "ExtGram")

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
end