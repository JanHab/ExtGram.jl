"""
    Example file for a 1D-1D Two-Stream instability test case.
"""

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LaTeXStrings

# Access arguments by index
M_vector = [4, 8, 12, 16, 20, 24, 28, 32, 40] #!parse.(Int, ARGS[1]) # Vector of M values to test
closure = "ExtGram" #!ARGS[2] # String # "Gram", "ExtGram" or "Grad"
T_end = 50 #!parse(Float64, ARGS[3])
base_tree_level = 8 #!parse(Int, ARGS[4]) # e.g. 10
polydeg = 1 #!parse(Int, ARGS[5]) # e.g. 1
ϵ = 0.01 #!parse(Float64, ARGS[6]) # 0.01
k = 0.5 #!parse(Float64, ARGS[7]) # 0.5

x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

Kn = 1.0 # doesn't matter, as zero-relaxation in the vlasov_poisson_source_term

ρ_0, v_0, θ_0 = 1.0, 0.0, 1.0 # initial conditions for Landau damping


for M in M_vector
    #= set up equations and initial condition =#
    equations = GramianMomentEquations1D(M, Kn, closure)
    initial_condition = InitialConditionsLandauDamping(
        ρ_0,
        ϵ,
        v_0,
        θ_0,
        k,
        equations
    )

    #= set up semidiscretization =#
    basis = LobattoLegendreBasis(polydeg)

    surface_flux = flux_lax_friedrichs
    volume_flux = flux_central
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

    solver = DGSEM(basis, surface_flux, volume_integral)

    mesh = TreeMesh(
        (domain[1],), (domain[2],), 
        initial_refinement_level=base_tree_level, n_cells_max=10_000, 
        periodicity=true
    )

    vlasov = VlasovPoissonCache(mesh, solver)   # buffers, FFT plans and field history for this run
    semi = SemidiscretizationHyperbolic(
        mesh, equations, 
        initial_condition, solver, 
        source_terms=VlasovPoissonSource(vlasov)
    )

    #= set up ODE =#
    tspan = (0.0, T_end)
    ode = semidiscretize(semi, tspan)

    callbacks, summary_callback = callbacksGramianMomentEquations(
        semi, tspan, basis; 
        cfl = 0.99,          # Maximum cfl number
        plot_interval = 20,  # plot every 20 steps
        time_interval = 500, # save at 500 time intervals
        name="VlasovPoisson/LandauDamping/moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name used for output
    )
    # Add Vlasov-Poisson callback
    callbacks = CallbackSet(callbacks, vlasov_poisson_callback(vlasov))

    #= solve =#
    sol = solve(
        ode, 
        CarpenterKennedy2N54(
            williamson_condition = false
        );
        dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
        ode_default_options()..., 
        save_everystep=true,
        callback = callbacks,
    );

    summary_callback()

    # Access the energy history
    E_L2_history = vlasov.E_L2
    time_callback = vlasov.times

    # store to csv file
    CSV.write(
        "out/VlasovPoisson/LandauDamping/energy_moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_p$(polydeg)_level$(base_tree_level)_x_lower$(x_lower)_x_upper$(x_upper).csv",
        Tables.columntable((
            time=time_callback, E_L2=E_L2_history, 
            E_L2_normalized=E_L2_history ./ E_L2_history[1], 
        ))
    )

    plot(
        time_callback, E_L2_history ./ E_L2_history[1], 
        xlabel = "Time", ylabel = L"$\|E\|_2 / \|E\|_2(t=0)$", 
        title = "Landau Damping: Electric Field Energy History",
        legend = false,
        lw = 2,
        grid = true,
        yscale = :log10,
    )
    savefig("Landau_M$(M).png")
end