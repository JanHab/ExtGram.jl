""" 
    Some helper functions for setting up 1D Gramian moment equations with Riemann initial conditions
"""

function setupGramianMomentEquations1DRiemann(
    M, Kn, closure, 
    f_left, f_right;
    base_tree_level = 8,
    surface_flux = flux_lax_friedrichs,
    volume_flux = flux_central,
    polydeg = 1,        # DG polynomial degree
    domain = (-2.0, 2.0),
    χ_set = "optimal",
    alpha_max = 0.5,
    alpha_min = 0.001,
    alpha_smooth = true,
    )
    """
        Setup ready to create semidiscretizations of the Gramian moment equations in 1D

    # Arguments
    - `M`: Number of moments
    - `Kn`: Knudsen number
    - `closure`: Closure function
    - `f_left`: Left distribution function for Riemann initial condition
    - `f_right`: Right distribution function for Riemann initial condition
    - `base_tree_level=8`: Base tree level for the mesh
    - `surface_flux=flux_lax_friedrichs`: Surface flux function
    - `volume_flux=flux_central`: Volume flux function
    - `polydeg=1`: Polynomial degree for DG basis
    - `domain=(-2.0, 2.0)`: Spatial domain
    - `χ_set="optimal"`: Set of χ values for closure - only relevant for Extended Gramian Closure with even M
    - `alpha_max=0.5`: Maximum shock capturing parameter
    - `alpha_min=0.001`: Minimum shock capturing parameter
    - `alpha_smooth=true`: Smooth shock capturing parameter
    """
    equations = GramianMomentEquations1D(M, Kn, closure; χ_set=χ_set)
    initial_condition = InitialConditionsShockTube(
        f_left, # Density, velocity, temperature
        f_right, # Shock in density, but not velocity, temperature initially
        equations
    )

    #= set up semidiscretization =#
    basis = LobattoLegendreBasis(polydeg)

    # shock capturing
    indicator_sc = IndicatorHennemannGassner(
        equations, basis,
        alpha_max = alpha_max,
        alpha_min = alpha_min,
        alpha_smooth = alpha_smooth, # smoothes with all neighboring indicators to remove numerical artifacts
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

    return basis, mesh, equations, initial_condition, solver, boundary_conditions
end


function callbacksGramianMomentEquations(
    semi, tspan, basis;
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    time_interval = 20, # save at 20 time intervals
    name="gram_solution",
)
    """
        Creates set of callbacks

    # Arguments
    - `semi`: Semidiscretization (Trixi.jl)
    - `tspan`: Time span of the simulation
    - `basis`: DG basis
    - `cfl=0.45`: Maximum CFL number
    - `plot_interval=20`: Plot every `plot_interval` steps
    - `time_interval=20`: Save solution `time_interval` often times
    - `name="gram_solution"`: Name of the output files (*.tsv)
    """
    alive_callback = AliveCallback(analysis_interval=100)
    summary_callback = SummaryCallback()
    stepsize_callback = StepsizeCallback(cfl=cfl)

    plot_callback = VisualizationCallback(
        semi;
        interval=plot_interval,
        solution_variables=cons2cons,
        plot_data_creator=PlotData1D,
        plot_creator=Trixi.show_plot
    )

    save_solution = SaveTriangulationCallback(
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
        save_solution,
        save_solution_prim,
    )

    return callbacks, summary_callback
end
