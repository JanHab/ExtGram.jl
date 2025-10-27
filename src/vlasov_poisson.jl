# Global storage for the electric field
# This will be updated during each RHS evaluation
mutable struct ElectricFieldStorage
    E::Vector{Float64}
    domain::Tuple{Float64, Float64}
    x_range::Vector{Float64}
    n::Int
    initialized::Bool
    ρ::Vector{Float64} # todo: remove this later
    MP1::Int
    E_L2::Vector{Float64}
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
    variables#::Matrix{Float64} # todo: remove this later
end

# Global instance
const ELECTRIC_FIELD = ElectricFieldStorage(Float64[], (0.0, 0.0), Float64[], 0, false, Float64[], 0, Float64[], Float64[],
    [Float64[]], # todo: remove this later
)

# Source term that solves Poisson globally and applies local source
function vlasov_poisson_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    # todo: hard-coded for the moment
    x_range = ELECTRIC_FIELD.x_range
    E_field = linear_interpolation(x_range, ELECTRIC_FIELD.E, extrapolation_bc = Interpolations.Line())
    
    # Evaluate electric field at position x
    E_local = E_field(x[1])

    # Apply source terms
    source = MVector{Mp1, Float64}(undef)
    source[1] = 0.0                    # No source for mass conservation

    # higher moments
    for k in 2:ELECTRIC_FIELD.MP1
        # Electric force affects momentum (2nd moment) and higher moments
        source[k] = (k-1) * E_local * u[k-1]
    end

    # return positive value of source, as it is on the lhs but with a minus -> positive source term
    # todo: Add possibility to combine with other source terms
    return source #!.+ relaxation_source(u, x, t, equations)
end

function solve_poisson_periodic_fft(ρ::AbstractVector{<:Real})
    n = length(ρ)
    Lx = ELECTRIC_FIELD.domain[end] - ELECTRIC_FIELD.domain[1]
    ρ̃  = ρ .- 1#! mean(ρ)  # neutralizing background
    ρk = fft(ρ̃)

    # build wavenumbers k consistent with FFT ordering
    # k = 0, 1, ..., floor(n/2), -ceil((n-1)/2), ..., -1
    k_int = [0:div(n,2); -div(n-1,2):-1]
    kx = (2π / Lx) .* k_int

    Ek = similar(ρk)
    Ek[1] = 0 # k = 0 mode
    for j in 2:n
        if kx[j] != 0.0
            # Directly solves for E in Fourier space
            Ek[j] = ρk[j] / (im * kx[j]) # ? Which sign is correct???
            # Alternative: solves for potential and then differentiate
            # ! Expression E is wrong here!!!
            # ? Ek[j] = ρk[j] / (kx[j]^2)
        else
            Ek[j] = 0
        end
    end

    E = real(ifft(Ek))
    return E
end

# Callback to solve Poisson equation globally at each timestep
function vlasov_poisson_callback(integrator)
    u = integrator.u
    t = integrator.t
    # todo: Use collect1dTreeArrays from converter_1d.jl
    semi = integrator.p
    u_ode = integrator.u

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
    ρ = variables[2:end-1, 1]  # First column is density # todo: why 2:end-1? What's wrong here? The first and last entry seem to be off, though.
    ELECTRIC_FIELD.variables = variables # todo: remove this later
    ELECTRIC_FIELD.ρ = vec(ρ)
    ELECTRIC_FIELD.x_range = vec(coordinates[2:end-1])

    # Solve Poisson equation globally
    E = solve_poisson_periodic_fft(ρ)
    # Store the electric field
    ELECTRIC_FIELD.E = copy(E)

    # Add energy vector
    # Energy = ||E(t,⋅)||_L2 = (∫ |E(t,x)|² dx)^(1/2)  (approximated via trapezoidal rule)
    # L2-norm of electric field
    E_L2 = trapz(ELECTRIC_FIELD.x_range, E.^2)^(1/2)
    # E_L2 = 1/2 * sum(E.^2) * (ELECTRIC_FIELD.x_range[2] - ELECTRIC_FIELD.x_range[1]) / n_cells
    push!(ELECTRIC_FIELD.E_L2, E_L2)
    push!(ELECTRIC_FIELD.times, t)  # Store the actual time

    ELECTRIC_FIELD.initialized = true
    
    return nothing
end

# Create the callback - triggers after each iteration
function vlasov_poisson_callback(;M, mesh, domain)
    # Reset the global storage to clear old data from previous runs
    empty!(ELECTRIC_FIELD.E)
    empty!(ELECTRIC_FIELD.E_L2)
    empty!(ELECTRIC_FIELD.times)
    empty!(ELECTRIC_FIELD.ρ)
    ELECTRIC_FIELD.n = 0
    ELECTRIC_FIELD.initialized = false
    ELECTRIC_FIELD.domain = domain
    
    # Set parameters for this run
    ELECTRIC_FIELD.MP1 = M+1
    # x_range = sort(mesh.tree.coordinates[1,1:mesh.tree.length])
    # x_range = (x_range[1:end-1] + x_range[2:end]) / 2  # cell centers
    # pushfirst!(x_range, domain[1])
    # push!(x_range, domain[2])
    # # x_range = [0, (x_coordinates[2:end] + x_coordinates[1:end-1])/2, 4.0*π]  # cell centers including boundaries
    # ELECTRIC_FIELD.x_range = x_range
    
    return DiscreteCallback(
        (u, t, integrator) -> true,  # Always trigger at every step
        vlasov_poisson_callback,
        save_positions=(true, true), # ?(false, false),
        initialize = (c, u, t, integrator) -> vlasov_poisson_callback(integrator)  # Call at initialization
    )
end

struct InitialConditionsCosine{N}
    ρ0::Float64
    ϵ::Float64
    v0::Float64
    θ0::Float64
    k::Float64

    function InitialConditionsCosine(ρ0::Float64, ϵ::Float64, v0::Float64, θ0::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(ρ0, ϵ, v0, θ0, k)
    end
end

function (ic::InitialConditionsCosine)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    ρx = ic.ρ0 * (1 + ic.ϵ * cos(ic.k * coords[1]))
    f = Maxwellian(ρx, ic.v0, ic.θ0)
    return convective_moments(f, Val(Mp1))
end
