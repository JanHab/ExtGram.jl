# plot first three moments
function plot_ρ_v_p(sol, M, x_lower, x_upper)
    u_final = sol.u[end]
    L = length(u_final)
    Nloc = L ÷ (M+1)
    Fmat = reshape(u_final, M+1, Nloc)

    ρ = Fmat[1, :]
    v = Fmat[2, :]
    ρΘ = Fmat[3, :]
    p = ρΘ

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
function TVD_space(sol, M::Int)
    L = length(sol.u[end]) # get number of local cells (from last time)
    Nloc = L ÷ (M+1)
    
    nt = length(sol.u)
    TV_t = zeros(M+1, nt)

    for n in 1:nt
        u_n = reshape(sol.u[n], M+1, Nloc)   # shape (M, nx)
        for m in 1:M+1
            TV_t[m, n] = sum(abs.(diff(u_n[m, :]))) # m-th moment
        end
    end

    p = plot()
    for m in 1:M+1
        plot!(
            p,
            sol.t, TV_t[m, :],
            marker=:o,
            label = "u^{($(m-1))}",
            xlabel="t",
            ylabel="TV in space",
            legend=:topright
        )
    end

    return p, TV_t
end

# read solutino from tvd file
function readsol(filename)
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
    x = Vector{Float64}
    u_solutions = Matrix{Float64}  # (vars × npts) matrix per dataset

    data_blocks = readsol(filename)              # Vector of matrices (npts × (1 + n_vars))
    final_block = data_blocks[end]               # last time level
    x = final_block[:, 1]                   # x coordinates (npts)
    y = final_block[:, 2:end]                    # moments (npts × n_vars)
    u_solutions = permutedims(y)            # shape: (n_vars × npts)

    return x, u_solutions
end
