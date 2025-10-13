# Setup ready to create semidiscretizations of the Gramian moment equations in 1D
function setupGramianMomentEquations(
    M, Kn, extended, 
    f_left, f_right;
    base_tree_level = 8,
    surface_flux = flux_lax_friedrichs,
    volume_flux = flux_central,
    polydeg = 1,        # DG polynomial degree
    domain = (-2.0, 2.0),
    )
    equations = GramianMomentEquations1D(M, Kn, extended)
    initial_condition = InitialConditionsShockTube(
        f_left, # Density, velocity, temperature
        f_right, # Shock in density, but not velocity, temperature initially
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

    return basis, mesh, equations, initial_condition, solver, boundary_conditions
end


# Creates set of callbacks
function callbacksGramianMomentEquations(
    semi, tspan, basis;
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
)
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

    # save_solution = SaveTriangulationCallback(
    #     time_interval=tspan[2]/20,
    #     save_initial_solution=true,
    #     file_format="tsv",
    #     append_solution=true,
    #     solution_variables = cons2cons,
    #     clear_out_dir=false,
    #     name=name,
    #     info="basis = $(Base.typename(typeof(basis)).wrapper)"
    # )

    callbacks = CallbackSet(
        alive_callback,
        stepsize_callback,
        plot_callback,
        # save_solution,
    )

    return callbacks, summary_callback
end
