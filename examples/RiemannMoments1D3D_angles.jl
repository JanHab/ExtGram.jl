using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra
using StaticArrays

# Access arguments by index
M = 4 #parse(Int, ARGS[1])
closure = "ExtGram" #ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = 1.0 #parse(Float64, ARGS[3])
# source_string = ARGS[4]
source = relaxation_source #zero_source #relaxation_source
T_end = 0.3 #parse(Float64, ARGS[5])
base_tree_level = 7 #!2 #parse(Int, ARGS[6]) # e.g. 8
polydeg = 1 #parse(Int, ARGS[7]) # e.g. 3
ρ_L = 7.0 #parse(Float64, ARGS[8]) # 7.0
v_L1 = 0.0 #!1.0 #!1.5 #parse(Float64, ARGS[9]) # 0.0
v_L2 = 0.0
v_L3 = 0.0
θ_L = 1.0 #parse(Float64, ARGS[10]) # 1.0
ρ_R = 1.0 #parse(Float64, ARGS[11]) # 1.0
v_R1 = 0.0 #!1.0 #parse(Float64, ARGS[12]) # 0.0
v_R2 = 0.0
v_R3 = 0.0
θ_R = 1.0 #parse(Float64, ARGS[13]) # 1.0


# Fixed settings
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)

angles1 = [ # maximizing angles
    (3.14159, 1.5708),
    (0.684719, 4.71239),
    (2.03444, 1.5708),
    (2.18628, 0.886077)
]
angles2 = [ # maximizing angles / 2
    (3.14159/2, 1.5708/2),
    (0.684719/2, 4.71239/2),
    (2.03444/2, 1.5708/2),
    (2.18628/2, 0.886077/2)
]
angles3 = [
    (2.4, 0.3),
    (0.3, 3.314159),
    (1.3, 1.5),
    (2.9, 0.7)
]
angles4 = [
    (2.8, 1.1),
    (0.9, 4.3),
    (2.5, 1.1),
    (2.7, 0.6)
]
angles5 = [
    (2.3, 3.7),
    (1.3, 2.7),
    (2.9, 0.1),
    (0.7, 0.1)
]

angle_pairs = [angles1, angles2, angles3, angles4, angles5]
# angle_pairs = [
#     [ # maximizing angles * 2
#     (2.3, 3.7),
#     (1.3, 2.7),
#     (2.9, 0.1),
#     (0.7, 0.1)#(2.18628, 0.886077)
# ]
# ]

for (anglepair, angles) in enumerate(angle_pairs)
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

    cfl = 0.99
    time_interval = 20
    name = "1D3D/1D3D_angles/gram_solution_1D3D_closure$(closure)_Kn$(Kn)_anglepairnumber$(anglepair)_Tend$(T_end)"

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

    # Plot first moments
    pl = plot();
    plot!(
        pl,
        title="Riemann Problem Moments 1D3D: M=$(M), closure=$(closure), Kn=$(Kn), source=$(source), T_end=$(T_end)",
        xlabel="x",
        ylabel="Primitive Variables",
    );
    x = range(x_lower, x_upper, length=Nloc);

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
        label="300",
        color=:purple,
    );

    plot!(
        pl,
        x, u[7, :], 
        label="U400",
        color=:brown,
    );

    display(pl)
end