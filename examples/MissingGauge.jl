using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Parameter
M_vector = [4, 6, 8, 10, 12]    # number of moments
extended = true # flag for extended gramian closure or "standard" closure
Kn = 1.0  # Knudsen number
T_end = 0.3
source = relaxation_source
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)

ρ_L = 7.0; v_L = 0.0; θ_L = 1.0
ρ_R = 1.0; v_R = 0.0; θ_R = 1.0

# χ value -> Gauge for χ = (n+1)/n with n = M/2
χ_vector = [-1.0, 0.0, 1.0]
# copy the convergence routine for each χ value and adapt the χ in the equations
for χ in χ_vector
    println("Running simulations with χ = $χ")
    for M in M_vector
        println("Running simulation with M = $M moments")
        # Setting up everything
        basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DRiemann(
            M, Kn, extended,
            Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
            Maxwellian(ρ_R, v_R, θ_R);
            domain = domain,
            χ_set = χ # set χ value explicitly
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
            name="Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)", # name of output files
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
        savefig(p1, "out/Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.pdf")

        # Store primitive variables in CSV file
        CSV.write(
            "out/Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.csv",
            Tables.columntable((x=x, rho=ρ, v=v, p=p))
        )

        # Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
        n_plots = 5
        p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
        display(p2)
        savefig(p2, "out/Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_λ_max.pdf")

        # Plot total variation in space over time
        p3, TV_t = TVD_space(sol, M)
        display(p3)
        savefig(p3, "out/Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_TV_space.pdf")
    end

    # Load the solutions and plot convergence
    # for M in M_vector
    plot_moments = 5
    conservative_plots = []; primitive_plots = []
    solution_conservative = []; solution_primitive = []
    for _ in 1:plot_moments
        push!(conservative_plots, plot())
        push!(primitive_plots, plot())
    end
    for M in M_vector
        solution_file = "out/Convergence/ChiValues/Chi$(Int(χ))/gram_solution_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R).tsv"

        x, conservative_moments = readfile(solution_file)
        primitive_moments = similar(conservative_moments)
        for i in 1:size(conservative_moments, 2)
            primitive_moments[:, i] = moment_cons2prim(conservative_moments[:, i])
        end
        push!(solution_conservative, conservative_moments)
        push!(solution_primitive, primitive_moments)

        for i in 1:plot_moments
            # conservative variables
            plot!(
                conservative_plots[i],
                x, conservative_moments[i, :], 
                label = "M = $M",
                xlabel = "x",
                ylabel = "u^{($i)}",
                legend = :topright
            )
            # primitive variables
            plot!(
                primitive_plots[i],
                x, primitive_moments[i, :], 
                label = "M = $M",
                xlabel = "x",
                ylabel = "w^{($i)}",
                legend = :topright
            )
        end
    end
    for i in 1:plot_moments
        display(conservative_plots[i])
        savefig(conservative_plots[i], "out/Convergence/ChiValues/Chi$(Int(χ))/conservative$(i)_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R).pdf")
        display(primitive_plots[i])
        savefig(primitive_plots[i], "out/Convergence/ChiValues/Chi$(Int(χ))/primitive$(i)_M$(M)_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R).pdf")
    end

    # relative L2-error
    errors_conservative = zeros(length(M_vector)-1, plot_moments)
    errors_primitive = zeros(length(M_vector)-1, plot_moments)
    for i in 1:(length(M_vector)-1)
        for j in 1:plot_moments
            errors_conservative[i, j] = norm(solution_conservative[i][j,:] .- solution_conservative[end][j,:]) / norm(solution_conservative[end][j,:])
            errors_primitive[i, j] = norm(solution_primitive[i][j,:] .- solution_primitive[end][j,:]) / norm(solution_primitive[end][j,:])
        end
    end

    l2_plot_conservative = plot(title="Relative L2-error conservative variables", xlabel="Number of moments M", ylabel="Relative L2-error", yscale=:log10, xticks=M_vector[1:end-1])
    l2_plot_primitive = plot(title="Relative L2-error primitive variables", xlabel="Number of moments M", ylabel="Relative L2-error", yscale=:log10, xticks=M_vector[1:end-1])
    for j in 1:plot_moments
        plot!(
            l2_plot_conservative,
            M_vector[1:end-1], errors_conservative[:, j],
            label = "u^{($j)}",
            marker = :o,
            linestyle = :solid
        )
        plot!(
            l2_plot_primitive,
            M_vector[1:end-1], errors_primitive[:, j],
            label = "w^{($j)}",
            marker = :o,
            linestyle = :solid
        )
    end
    display(l2_plot_conservative)
    savefig(l2_plot_conservative, "out/Convergence/ChiValues/Chi$(Int(χ))/L2_error_conservative__Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R).pdf")
    display(l2_plot_primitive)
    savefig(l2_plot_primitive, "out/Convergence/ChiValues/Chi$(Int(χ))/L2_error_primitive__Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R).pdf")

end