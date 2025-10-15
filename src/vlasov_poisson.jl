# Global storage for the electric field
# This will be updated during each RHS evaluation
mutable struct ElectricFieldStorage
    E::Vector{Float64}
    x_coords::Vector{Float64}
    domain::Tuple{Float64, Float64}
    n::Int
    initialized::Bool
    counter::Int    # todo: remove this later
    ρ::Vector{Float64} # todo: remove this later
    MP1::Int
    Energy::Vector{Float64}  # todo: remove this later
end

# Global instance
# todo: make the domain dynamic
const ELECTRIC_FIELD = ElectricFieldStorage(Float64[], Float64[], (-0.0, 0.0), 0, false, 0, Float64[], 0, Float64[])

# Enhanced Poisson solver with proper boundary conditions
function solve_poisson_global(ρ, domain)
    n = length(ρ)
    x_min, x_max = domain
    Δx = (x_max - x_min) / (n - 1)
    
    # Solve ∂E/∂x = ρ by integration
    E = zeros(n)
    
    # Apply boundary condition: E = 0 at left boundary
    # E[1] = 0.0
    
    # # Integrate using trapezoidal rule
    # for i in 2:n
    #     E[i] = E[i-1] + 0.5 * (ρ[i] + ρ[i-1]) * Δx
    # end

    # Alternatively, solve ∂²ϕ/∂x² = -ρ with ϕ=0 at left boundary and ∂ϕ/∂x=0 at right boundary
    # Method is: Jacobi iteration
    ϕ = zeros(n)
    n_iter = 10000
    for iter in 1:n_iter
        for i in 2:n-1
            ϕ[i] = 0.5 * (ϕ[i-1] + ϕ[i+1] + ρ[i] * Δx^2)
        end
        # Boundary conditions for potential
        ϕ[1] = 0.0          # Dirichlet at left boundary
        ϕ[n] = 0.0 #!       # Dirichlet at right boundary ϕ[n-1]       # Neumann at right boundary (∂ϕ/∂x = 0)
    end

    # Compute electric field E = -∂ϕ/∂x central differences
    for i in 2:n-1
        E[i] = -(ϕ[i+1] - ϕ[i-1]) / (2 * Δx)
    end
    # Forward difference at left boundary
    E[1] = -(ϕ[2] - ϕ[1]) / Δx
    # Backward difference at right boundary
    E[n] = -(ϕ[n] - ϕ[n-1]) / Δx

    return E
end

# Source term that solves Poisson globally and applies local source
function vlasov_poisson_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    
    # todo: hard-coded for the moment
    x_min, x_max = ELECTRIC_FIELD.domain
    x_range = range(x_min, x_max, length=ELECTRIC_FIELD.n)
    E_field = interpolate((x_range,), ELECTRIC_FIELD.E, Gridded(Linear()))
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

# Callback to solve Poisson equation globally at each timestep
function vlasov_poisson_callback(integrator)
    u = integrator.u
    t = integrator.t
    
    # Extract density from the full solution vector
    MP1 = ELECTRIC_FIELD.MP1
    n_cells = length(u) ÷ MP1
    
    # Reshape to extract density
    U_matrix = reshape(u, MP1, n_cells)
    ρ = U_matrix[1, :]  # First row is density

    # Solve Poisson equation globally
    # E = solve_poisson_global(ρ .- 1, ELECTRIC_FIELD.domain) # todo: change back
    E = solve_poisson_global(ρ, ELECTRIC_FIELD.domain) # todo: change back

    # Store the electric field
    ELECTRIC_FIELD.E = copy(E)
    # Add energy vector
    # Energy = ||E(t,⋅)||_L2 = (∫ |E(t,x)|² dx)^(1/2)  (approximated via trapezoidal rule)
    Energy = (sum(E.^2) * (ELECTRIC_FIELD.domain[2] - ELECTRIC_FIELD.domain[1]) / n_cells)^(1/2)
    push!(ELECTRIC_FIELD.Energy, Energy)  # todo: remove this later
    # todo: remove below, just for debuggin purposes
    ELECTRIC_FIELD.ρ = ρ
    ELECTRIC_FIELD.x_coords = range(ELECTRIC_FIELD.domain[1], ELECTRIC_FIELD.domain[2], length=n_cells)
    ELECTRIC_FIELD.counter += 1    # todo: remove this later

    ELECTRIC_FIELD.n = n_cells
    ELECTRIC_FIELD.initialized = true
    
    return nothing
end

# Create the callback - triggers after each iteration
function vlasov_poisson_callback(;M, domain)
    # Option 1: DiscreteCallback that triggers at every accepted step
    ELECTRIC_FIELD.MP1 = M+1
    ELECTRIC_FIELD.domain = domain
    return DiscreteCallback(
        (u, t, integrator) -> true,  # Always trigger at every step
        vlasov_poisson_callback,
        save_positions=(false, false),
        initialize = (c, u, t, integrator) -> vlasov_poisson_callback(integrator)  # Call at initialization
    )
end
