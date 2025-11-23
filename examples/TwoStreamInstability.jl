using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LaTeXStrings, FastGaussQuadrature

# Parameter
M_vector = [4] #![4, 12, 24]    # number of moments
closure = "Gram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0 # ! doesn't matter, as zero-relaxation in the vlasov_poisson_source_term # Knudsen number
T_end = 25#!50.0
source = vlasov_poisson_source
x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

base_tree_level = 8 # ! 6#!5
polydeg = 3

for M in M_vector
    equations = GramianMomentEquations1D(M, Kn, closure)
        ϵ = 0.01
        k = 0.5
        initial_condition = InitialConditionsTwoStream(
            ϵ, k,
        equations
    )

    #= set up semidiscretization =#
    basis = LobattoLegendreBasis(polydeg)

    surface_flux = flux_lax_friedrichs #flux_central #!flux_lax_friedrichs
    volume_flux = flux_central #!flux_central
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

    solver = DGSEM(basis, surface_flux, volume_integral)

    mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=true)

    semi = SemidiscretizationHyperbolic(
        mesh, equations, 
        initial_condition, solver, 
        # boundary_conditions=boundary_conditions, 
        source_terms=source
    )

    #= set up ODE =#
    tspan = (0.0, T_end)
    ode = semidiscretize(semi, tspan)

    # callbacks, summary_callback = callbacksGramianMomentEquations(
    #     semi, tspan, basis; 
    #     cfl = 0.99,          # Maximum cfl number
    #     plot_interval = 20,  # plot every 20 steps
    #     time_interval = 500, # save at 500 time intervals
    #     # * name="VlasovPoisson/VlasovPoisson_moments_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level)", # name used for output
    #     name="VlasovPoisson/gram_moments", # name used for output
    # )
    alive_callback = AliveCallback(analysis_interval=100)
    summary_callback = SummaryCallback()
    cfl = 0.99
    stepsize_callback = StepsizeCallback(cfl=cfl)

    plot_callback = VisualizationCallback(
        semi;
        interval=20,
        solution_variables=cons2cons,
        plot_data_creator=PlotData1D,
        plot_creator=Trixi.show_plot
    )

    save_solution = SaveTriangulationCallback(
        time_interval=tspan[2]/20,
        save_initial_solution=true,
        file_format="tsv",
        append_solution=true,
        solution_variables = cons2cons,
        clear_out_dir=false,
        name=name="VlasovPoisson/TwoStreamInstability/moments_closure$(closure)_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level)_x_lower$(x_lower)_x_upper$(x_upper)", # name used for output
        info="basis = $(Base.typename(typeof(basis)).wrapper)"
    )

    callbacks = CallbackSet(
        alive_callback,
        stepsize_callback,
        plot_callback,
        save_solution,
    )
    # Add Vlasov-Poisson callback
    callbacks = CallbackSet(callbacks, vlasov_poisson_callback(;M, mesh, domain))

    #= solve =#
    # ? may need a fixed time-step to get the instabilities through.
    # ? Currently stopping when reaching it
    sol = solve(
        ode, 
        CarpenterKennedy2N54(
            williamson_condition = false
        );
        dt = 1.0,
        # Euler(); #Euler();# Vern6(); #Euler();
        # dt = 1/50, # solve needs some value here but it will be overwritten by the stepsize_callback
        ode_default_options()..., 
        save_everystep=true,
        callback = callbacks,
    );

    summary_callback()

    # Access the energy history
    E_L2_history = HyQMOM.ELECTRIC_FIELD.E_L2;
    time_callback = HyQMOM.ELECTRIC_FIELD.times;  # Use actual times from callback

    # You can then plot it or analyze it
    γ = -0.1533; # theoretical decay rate for k=1/2
    γt = exp.(γ .* time_callback);
    plot(
        time_callback, E_L2_history ./ E_L2_history[1],
        xlabel="Time", 
        label="HyQMOM M=$M (from callback)",
        ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
        yaxis=:log,
        legend=:bottomleft
    )
    savefig("out/VlasovPoisson/TwoStreamInstability/energy_from_callback.pdf")
    # store to csv file
    CSV.write(
        "out/VlasovPoisson/TwoStreamInstability/energy_moments_closure$(closure)_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level).csv",
        Tables.columntable((
            time=time_callback, E_L2=E_L2_history, 
            E_L2_normalized=E_L2_history ./ E_L2_history[1], 
            theoretical_decay=γt
        ))
    )
end


# # Plot IC
# x_vals = range(x_lower, x_upper, length=100)
# c_vals = range(-5.0, 5.0, length=100)

# f_ic(x, c) = 1/sqrt(2*π) * exp(-c^2 / 2) * c.^2 .* (1 + ϵ * cos(k * x))

# heatmap(
#     x_vals, c_vals, (x,c)->f_ic(x,c),
#     xlabel="x", ylabel="c",
#     title="Initial Condition f(x,c)",
#     colorbar_title="f(x,c)"
# )
# contourf(
#     x_vals, c_vals, (x,c)->f_ic(x,c),
# )
# savefig("out/VlasovPoisson/TwoStreamInstability/initial_condition_two_stream.pdf")

# function collect1DTreeArrays_local(semi, u_ode, solution_variables)
#     mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
#     tree = mesh.tree

#     n_active_cells = Trixi.count_leaf_cells(tree)
#     active_cell_ids = Trixi.leaf_cells(tree)

#     basis_nodes = solver.basis.nodes
#     n_vertices_per_cell = length(basis_nodes) # Trixi.nnodes(solver)
#     # total vertices amount (overlapping vertices combined) = #cell_inner_vertices + #cell_outer_vertices(=n_cells+1)
#     n_vertices = (n_vertices_per_cell-1)*n_active_cells+1

#     n_vars = Trixi.nvariables(equations)
#     u = Trixi.wrap_array_native(u_ode, mesh, equations, solver, cache)

#     # i_active_cell -> vertex ids
#     # in 1D a simple relation suffices: active_cells_to_vertices(i_cell) = i_cell*(n_vertices_per_cell-1)-3 .+ collect(1:n_vertices_per_cell)
#     # use a matrix anyway for ease of use and compatibility (no need to pass along n_active_cells for iteration over cells)
#     active_cells_to_vertices = zeros(Int, n_active_cells, n_vertices_per_cell)

#     # auxiliary arrays for the solution output
#     coordinates = zeros(n_vertices)
#     variables = zeros(n_vertices, n_vars)


#     # go through all active cells from left to right and compute the nodal coordinates
#     for i_cell=1:n_active_cells
#         i_cell_global = active_cell_ids[i_cell]
#         vertices = (i_cell-1)*(n_vertices_per_cell-1) .+ collect(1:n_vertices_per_cell)
#         active_cells_to_vertices[i_cell, :] = vertices
#         coordinates[vertices] = basis_nodes*0.5*Trixi.length_at_cell(tree, i_cell_global) .+ Trixi.cell_coordinates(tree, i_cell_global)[1]
#     end

#     # collect all nodal variable values
#     for i_var=1:n_vars
#         data = vec(u[i_var, .., :])
#         for i_cell=1:n_active_cells
#             index = 1 + (i_cell-1)*n_vertices_per_cell # DG -> use a non-overlapping vertex enumeration | alternatively increment indices
#             nodal_values = data[index:index+n_vertices_per_cell-1]
#             nodal_values[[1,end]] *= 0.5 # average the cell boundary values
#             variables[active_cells_to_vertices[i_cell, :], i_var] += nodal_values
#         end
#     end
#     # note: for a periodic boundary, the leftmost and rightmost vertices are NOT averaged
#     # undo the 0.5 factor of the leftmost and rightmost vertices, as they are not averaged
#     variables[[1,end], :] *= 2

#     # convert each nodes variables if needed
#     if solution_variables != cons2cons
#         n_out_vars = size(solution_variables(zeros(n_vars), equations))[1]
#         if n_vars == n_out_vars # in place for same output size as input size
#             for i=1:n_vertices; variables[i,:] = solution_variables(variables[i,:], equations); end
#         else
#             aux = zeros(n_vertices, n_out_vars)
#             for i=1:n_vertices; aux[i,:] = solution_variables(variables[i,:], equations); end
#             variables = aux
#         end
#     end

#     # create the connectivity for linear FE elements, given by [i, i+1] in each column, since vertices/coordinates are sorted
#     # note: just output linear elements in every case for maximum compatibility/support while the outputs remain consistent
#     connectivity = collect([1:n_vertices-1 2:n_vertices]') # n_elements = n_vertices-1

#     return connectivity, coordinates, variables
# end


# for i in 1:800:length(sol.t)
#     println("Time step $(i)/$(length(sol.t)): t=$(sol.t[i])")
#     _, coords, vars = collect1DTreeArrays_local(semi, sol.u[i], cons2cons)

#     x = coords[2:end-1]  # exclude ghost cells
#     c_vec = range(-3.0, 3.0, length=100)
#     ρ = vars[2:end-1, 1]
#     # plot(x, ρ,
#     #     xlabel="x", ylabel="ρ(x)",
#     #     title="Density ρ(x) at t=$(sol.t[i])",
#     # )
#     v = vars[2:end-1, 2] ./ ρ
#     θ = (vars[2:end-1, 3] ./ ρ) .- v.^2

#     f_grad = []
#     for i in 1:length(coords)-2
#         α = HyQMOM.solve_alpha(vars[i+1, :], ρ[i], v[i], θ[i])
#         nquad = 2*(M+2)
#         xF, wF = gausshermite(nquad)
#         s = sqrt(2 * θ[i])
#         c = v[i] .+ s .* xF
#         pref = ρ[i] / sqrt(2 * π * θ[i])
#         f(c) = pref * exp(-0.5 * ((c - v[i]) / s)^2 + sum(α[j] * c^(j-1) for j in 1:length(α)))
#         f_grad_ = []
#         for cj in c_vec
#             push!(f_grad_, f(cj))
#         end
#         push!(f_grad, f_grad_)
#     end
#     f_grad = hcat(f_grad...)'

#     heatmap(
#         x, c_vec, (x,c)->f_grad[findfirst(isequal(x), x), findfirst(isequal(c), c_vec)],
#         xlabel="x", ylabel="c",
#         title="Distribution function f(x,c) at t=$(sol.t[i]) reconstructed from Grad closure",
#         colorbar_title="f(x,c)",
#         # clims=(0.0, 1.0)
#     )
#     display(current())
#     # savefig("f_distribution_grad_t$(round(sol.t[i], digits=2)).pdf")
# end



# # Plot ρ, φ, E
# _, coords, vars = collect1DTreeArrays_local(semi, sol.u[end], cons2cons)

# x = coords[2:end-1]  # exclude ghost cells
# c_vec = range(-3.0, 3.0, length=100)
# ρ = vars[2:end-1, 1]
# E = HyQMOM.ELECTRIC_FIELD.E

# using Statistics
# plt = plot(
#     xlabel="x",
# )
# plot!(
#     plt,
#     x, ρ .- mean(ρ),
#     label="ρ - ρ̄"
# )
# plot!(
#     plt,
#     x, E,
#     label="E"
# )

# init_ρ, init_v, init_θ = [], [], []
# using FastGaussQuadrature
# Mp1 = M+1
# for x_ in x
#     max(c) = 1/sqrt(2*π) * exp(-c^2 / 2) * c^2 .* (1 + ϵ * cos(k * x_))
#     ξ, w = gausshermite(Mp1+1) # +1 for good measure, should not be necessary
#     C = sqrt(2*1.0) .* ξ; c = C .+ 0.0
#     fw = max.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*1.0)
#     ic = SVector{Mp1,Float64}(ntuple(n->sum(c .^(n-1) .* fw), Mp1))
#     push!(init_ρ, ic[1])
#     push!(init_v, ic[2])
#     push!(init_θ, ic[3])
# end

# plt = plot(layout=3)
# plot!(
#     plt, subplot=1,
#     x, init_ρ,
#     label="ρ",
#     ylims=(0.98,1.02)
# )
# plot!(
#     plt, subplot=2,
#     x, init_v,
#     label="v",
#     ylims=(-0.01,0.01)
# )
# plot!(
#     plt, subplot=3,
#     x, init_θ,
#     label="θ",
#     ylims=(2.95,3.05)
# )