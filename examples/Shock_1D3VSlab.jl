
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays

M = parse(Int, ARGS[1])
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = parse(Float64, ARGS[3])
source_string = ARGS[4]
source = source_string == "relaxation_source" ? relaxation_source : zero_source 
angle_name = ARGS[5] # "Max", "Arc", "Det", "Bad"
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
domain = (-2.0, 2.0)
base_tree_level = 10
T_end = 0.3 
cfl = 0.45
time_interval = 10 # ? Increase later? 10
slab_geometry = true # ? false


theta_Max = [
    [3.14159, 0.616578, 1.10715, 2.18628, ], # M=4
    [3.14155, 2.80135, 0.857072, 0.684719, 1.1832, 2.0944], # M=6
    [3.14159, 0.390413, 0.0000975558, 2.57765, 3.14159, 0.785398, 1.91063, 1.10711, 1.10715] # M=8
]
phi_Max = [
    [1.5708, 4.39298, 4.71239, 0.886077], # M=4
    [1.5708, 2.02967, 4.71239, 5.27633, 4.71239, 2.28452], # M=6
    [1.5708, 5.01795, 3.98292, 1.07991, 2.52616, 5.32787, 1.5708, 5.44221, 3.98266], # M=8
]

theta_Max_odd = [
    [3.14156, 0.39418, 0.955317, 0.785398, 9.07658e-316, 1.5708], # M=5
    [3.14159, 0.357662, 3.14159, 2.52611, 2.0944, 0.886077, 6.95333e-310, 1.5708, 1.57027], # M=7
    [3.14159, 0.340691, 3.14159, 0.523599, 0., 2.42787, 2.03444, 0.955287, 2.18628, 3.25075e-315, 1.5708, 1.5708], # M=9
]

phi_Max_odd = [
    [1.5708, 2.05657, 1.5708, 0.955317, 9.07658e-316, 0.785398], # M=5
    [1.5708, 1.17981, 5.49779, 4.18879, 1.5708, 2.22985, 1.5e-323, 2.35624, 0.785253], # M=7
    [1.5708, 4.39155, 0.886077, 1.10715, 5.59824, 5.29203, 4.71239, 4.02762, 2.25552, 3.25082e-315, 0.477658, 0.785398], # M=9
]

theta_Arc = [
    [5.85555e-10, 3.14159, 3.14159, 1.57601], # M=4
    [0., 3.14159, 3.14159, 0.593484, 3.14159, 1.57352], # M=6
    [0., 0., 0., 0.386383, 3.14159, 0.785703, 0., 1.56657, 1.56918], # M=8
]

phi_Arc = [
    [4.71239, 5.30183, 0.00368397, 2.35619], # M=4
    [1.5708, 1.18628, 3.78509, 1.06088, 0.0024264, 3.92699], # M=6
    [1.5708, 4.99513, 5.38713, 4.35201, 3.60512, 2.18642, 0.00181229, 5.94573, 2.35619] # M=8
]

theta_Det = [
    [1.08562, 0.969517, 0.0841159, 0.140558], # M=4
    [0.791178, 2.71441, 1.19527, 0.7375, 2.69368, 2.97235], # M=6
    [0.787033, 2.64384, 0.0612564, 2.99817, 0.313342, 0.516051, 2.44407, 0.0226567, 0.99256], # M=8
]

phi_Det = [
    [1.57834, 2.11307, 4.62546, 4.0341], # M=4
    [1.49439, 5.46196, 0.779535, 0.457368, 4.84181, 1.53241], # M=6
    [5.32313, 4.40688, 5.47647, 5.13619, 0.599894, 5.59516, 2.9642, 1.72031, 2.16574], # M=8
]

theta_Bad = [ # M=4
    0.0, 1.0, 2.0, 3.0
]
phi_Bad = [ # M=4
    π*1/4, π*2/4, π*3/4, π*4/4
]

rho_L_vec = [7.0, 15.0, 30.0]
if closure == "Grad"
    rho_L_vec = [7.0]
end
for rho_L in rho_L_vec
    if angle_name == "Max"
        if isodd(M)
            theta = theta_Max_odd[M ÷ 2 - 1]
            phi = phi_Max_odd[M ÷ 2 - 1]
        else
            theta = theta_Max[M ÷ 2 - 1]
            phi = phi_Max[M ÷ 2 - 1]
        end
    elseif angle_name == "Arc"
        if isodd(M)
            error("Arc angles are not defined for odd M.")
        end
        theta = theta_Arc[M ÷ 2 - 1]
        phi = phi_Arc[M ÷ 2 - 1]
    elseif angle_name == "Det"
        if isodd(M)
            error("Det angles are not defined for odd M.")
        end
        theta = theta_Det[M ÷ 2 - 1]
        phi = phi_Det[M ÷ 2 - 1]
    elseif (angle_name == "Bad" && M == 4)
        theta = theta_Bad
        phi = phi_Bad
    else
        error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", \"Det\".")
    end

    equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=slab_geometry, theta=theta, phi=phi)

    f_left = Maxwellian1D3V(rho_L, (0.0, 0.0, 0.0), 1.0)
    f_right = Maxwellian1D3V(1.0, (0.0, 0.0, 0.0), 1.0)

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

    name = "1D3V/1D3V_angles/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_anglepair$(angle_name)_rhoL$(rho_L)"
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