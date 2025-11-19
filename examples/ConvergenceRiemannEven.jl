using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Parameter
M = 4    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.3
source = relaxation_source
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)
polydeg = 2
collect_tree_levels = LinRange{Int}(2, 10, 9) # levels to be tested

ρ_L = 7.0; v_L = 0.0; θ_L = 1.0
ρ_R = 1.0; v_R = 0.0; θ_R = 1.0

for tree_level in collect_tree_levels
    println("Running simulation with base refinement level = $tree_level")
    # Setting up everything
    basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DRiemann(
        M, Kn, closure,
        Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
        Maxwellian(ρ_R, v_R, θ_R);
        domain = domain,
        polydeg = polydeg,
        base_tree_level = tree_level,
    )

    semi = SemidiscretizationHyperbolic(
        mesh, equations, 
        initial_condition, solver, 
        boundary_conditions=boundary_conditions, 
        source_terms=source
    )

    #= set up ODE =#
    tspan = (0.0, T_end)
    ode = semidiscretize(semi, tspan)

    callbacks, summary_callback = callbacksGramianMomentEquations(
        semi, tspan, basis; 
        cfl = 0.45,          # Maximum cfl number
        plot_interval = 20,  # plot every 20 steps
        name="Convergence/ConvergenceEvenOrder/gram_solution_treelevel$(tree_level)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)", # name of output files
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
        saveat = range(tspan[1], tspan[2], length=100)
    );

    summary_callback()

    # Post Processing
    x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)
    display(p1)
    savefig(p1, "out/Convergence/ConvergenceEvenOrder/gram_solution_treelevel$(tree_level)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.pdf")

    # Store primitive variables in CSV file
    CSV.write(
        "out/Convergence/ConvergenceEvenOrder/gram_solution_treelevel$(tree_level)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.csv",
        Tables.columntable((x=x, rho=ρ, v=v, p=p))
    )
end
