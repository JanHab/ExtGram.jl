"""
    This file provides analysis scripts for the 1D1D moment equations
    
    Including: 
        - plot_ρ_v_p: calculates and plots density, velocity and pressure from the solution
        - plot_λ_max: calculates and plots the maximum eigenvalue of the flux Jacobian over time
        - conservation: calculates and plots the total variation in space of density, momentum and energy over time
        - readsol: reads solution data from a .tvd file (output)
        - readfile: reads the final time level solution from a .tvd file
"""

function plot_ρ_v_p(sol, M, x_lower, x_upper)
    """
    Calculates discretization, density, velocity and pressure

    Careful: Only works for polydeg=1 (linear basis functions)!

    # Arguments:
    - `sol`: solution object from Trixi.jl
    - `M`: maximum moment degree
    - `x_lower`: lower bound of the spatial domain
    - `x_upper`: upper bound of the spatial domain

    # Returns:
    - `x`: spatial discretization points
    - `ρ`: density values at discretization points
    - `v`: velocity values at discretization points
    - `p`: pressure values at discretization points
    - `p_plot`: Plots.jl plot object containing the plots of ρ, v, and p
    """
    u_final = sol.u[end]
    L = length(u_final)
    Nloc = L ÷ (M+1)
    Fmat = reshape(u_final, M+1, Nloc)

    ρ = Fmat[1, :]
    v = Fmat[2, :] ./ ρ
    ρΘplusρv2 = Fmat[3, :]
    ρθ = ρΘplusρv2 - ρ .* v .^ 2
    p = ρθ

    x = range(x_lower, x_upper, length=Nloc)

    p_plot = Plots.plot()
    plot!(
        p_plot,
        x, ρ,
        color=:blue,
        label="Density",
        legend=:topleft
    )
    plot!(
        twinx(p_plot),
        x, v, 
        color=:red,
        label="Velocity",
        legend=:topright,
    )
    plot!(
        p_plot,
        x, p,
        color=:green,
        label="Pressure",
    )

    return x, ρ, v, p, p_plot
end

# Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
function plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
    """
    Plots the maximum eigenvalue (wave-speed) of the flux Jacobian over time.

    Careful: Only works for polydeg=1 (linear basis functions)!

    # Arguments:
    - `semi`: semi-discretization object from Trixi.jl
    - `sol`: solution object from Trixi.jl
    - `M`: maximum moment degree
    - `n_plots`: number of time levels to plot
    - `x_lower`: lower bound of the spatial domain
    - `x_upper`: upper bound of the spatial domain

    # Returns:
    - `p`: Plots.jl plot object containing the maximum eigenvalue plots
    """
    L = length(sol.u[end]) # get number of local cells (from last time)
    Nloc = L ÷ (M+1)

    # x-domain (assumes mesh from -2 to 2)
    # todo: make this generic
    x = range(x_lower, x_upper, length=Nloc)

    # Plotting
    p = Plots.plot()
    nt = length(sol.u)
    n_use = min(n_plots, nt)
    indices = [round(Int, 1 + (k-1)*(nt-1)/(n_use-1)) for k in 1:n_use]  # equally spaced, includes 1 and nt
    for j in indices
        u = sol.u[j]
        Fmat = reshape(u, M+1, Nloc)
        λ_max = similar(Fmat, Nloc)
        for i in 1:Nloc
            λ_max[i] = maximum(abs.(real.(eigen(flux_jacobian(Fmat[:, i], semi.equations)).values)))
        end
        plot!(
            p,
            x, λ_max,
            label = "t = $(round(sol.t[j], digits=3))",
            legend = :topleft
        )
        plot!(p, xlabel="x", ylabel="λ_max")
    end
    return p
end

# Plot total variation in space of density over time
function conservation(sol, M::Int, semi)
    """
    Plots the total variation in space of density, momentum and energy over time.

    # Arguments:
    - `sol`: solution object from Trixi.jl
    - `M::Int`: maximum moment degree
    - `semi`: semi-discretization object from Trixi.jl

    # Returns:
    - `p`: Plots.jl plot object containing the total variation plots
    - `mass`: array of total mass over time
    - `momentum`: array of total momentum over time
    - `energy`: array of total energy over time
    """
    solution_variables = cons2cons

    mass, momentum, energy = [], [], []
    for u_ode in sol.u
        mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
        tree = mesh.tree

        n_active_cells = Trixi.count_leaf_cells(tree)
        active_cell_ids = Trixi.leaf_cells(tree)

        basis_nodes = solver.basis.nodes
        n_vertices_per_cell = length(basis_nodes) # Trixi.nnodes(solver)
        # total vertices amount (overlapping vertices combined) = #cell_inner_vertices + #cell_outer_vertices(=n_cells+1)
        n_vertices = (n_vertices_per_cell-1)*n_active_cells+1

        n_vars = Trixi.nvariables(equations)
        u = Trixi.wrap_array_native(u_ode, mesh, equations, solver, cache)

        # i_active_cell -> vertex ids
        # in 1D a simple relation suffices: active_cells_to_vertices(i_cell) = i_cell*(n_vertices_per_cell-1)-3 .+ collect(1:n_vertices_per_cell)
        # use a matrix anyway for ease of use and compatibility (no need to pass along n_active_cells for iteration over cells)
        active_cells_to_vertices = zeros(Int, n_active_cells, n_vertices_per_cell)

        # auxiliary arrays for the solution output
        coordinates = zeros(n_vertices)
        variables = zeros(n_vertices, n_vars)


        # go through all active cells from left to right and compute the nodal coordinates
        for i_cell=1:n_active_cells
            i_cell_global = active_cell_ids[i_cell]
            vertices = (i_cell-1)*(n_vertices_per_cell-1) .+ collect(1:n_vertices_per_cell)
            active_cells_to_vertices[i_cell, :] = vertices
            coordinates[vertices] = basis_nodes*0.5*Trixi.length_at_cell(tree, i_cell_global) .+ Trixi.cell_coordinates(tree, i_cell_global)[1]
        end

        # collect all nodal variable values
        for i_var=1:n_vars
            data = vec(u[i_var, .., :])
            for i_cell=1:n_active_cells
                index = 1 + (i_cell-1)*n_vertices_per_cell # DG -> use a non-overlapping vertex enumeration | alternatively increment indices
                nodal_values = data[index:index+n_vertices_per_cell-1]
                nodal_values[[1,end]] *= 0.5 # average the cell boundary values
                variables[active_cells_to_vertices[i_cell, :], i_var] += nodal_values
            end
        end
        # note: for a periodic boundary, the leftmost and rightmost vertices are NOT averaged
        # undo the 0.5 factor of the leftmost and rightmost vertices, as they are not averaged
        variables[[1,end], :] *= 2

        # convert each nodes variables if needed
        if solution_variables != cons2cons
            n_out_vars = size(solution_variables(zeros(n_vars), equations))[1]
            if n_vars == n_out_vars # in place for same output size as input size
                for i=1:n_vertices; variables[i,:] = solution_variables(variables[i,:], equations); end
            else
                aux = zeros(n_vertices, n_out_vars)
                for i=1:n_vertices; aux[i,:] = solution_variables(variables[i,:], equations); end
                variables = aux
            end
        end

        ρ = variables[2:end-1, 1]
        ρv = variables[2:end-1, 2]
        ρΘplusρv2 = variables[2:end-1, 3]

        push!(mass, trapz(coordinates[2:end-1], ρ))
        push!(momentum, trapz(coordinates[2:end-1], ρv))
        push!(energy, trapz(coordinates[2:end-1], ρΘplusρv2))
    end

    # nt = length(sol.u)
    # TV_t = zeros(3, nt)

    p = plot(
        layout = (3, 1),
        label= "t",
        legend = false,
    )
    for (i, (property, label)) in enumerate(zip([mass, momentum, energy], ["Density", "Momentum", "Energy"]))
        plot!(
            p,
            # sol.t, property,
            sol.t[2:end], diff(property),
            marker=:o,
            label = label,
            legend=:topright,
            subplot = i
        )
    end
    # plot!(p, yscale=:log10)

    return p, mass, momentum, energy
end

# read solution from tvd file
function readsol(filename)
    """
    Reads solution data from a .tvd file (output).

    # Arguments:
    - `filename`: path to the .tvd file

    # Returns:
    - `blocks`: Vector of matrices, each matrix corresponds to a time level
    """
    f = open(filename)
    lines = readlines(f)
    n_vars = length(split(lines[3]))-1 # counts coordinates x as var

    # split into time-level blocks
    blocks = []
    i = 1
    while i < length(lines)
        l = lines[i]
        if startswith(l, "# timestep")# && endswith(l, "t=$T_end")
        j = findfirst("points=", l).stop
        n_points = parse(Int, l[j+1:j+findfirst(",", l[j:end]).start-2])

        data = zeros(Float64, n_points, n_vars)
        for k=1:n_points
            data[k, :] .= parse.(Float64, split(lines[i+k]))
        end
        push!(blocks, data)
        i += n_points
        else i += 1 end
    end

    return blocks
end

function readfile(filename)
    """
    Reads the final time level solution from a .tvd file.

    # Arguments:
    - `filename`: path to the .tvd file

    # Returns:
    - `x`: Vector of spatial coordinates
    - `u_solutions`: Matrix of solution variables (variables × npts)
    """
    x = Vector{Float64}
    u_solutions = Matrix{Float64}  # (vars × npts) matrix per dataset

    data_blocks = readsol(filename)              # Vector of matrices (npts × (1 + n_vars))
    final_block = data_blocks[end]               # last time level
    x = final_block[:, 1]                   # x coordinates (npts)
    y = final_block[:, 2:end]                    # moments (npts × n_vars)
    u_solutions = permutedims(y)            # shape: (n_vars × npts)

    return x, u_solutions
end
