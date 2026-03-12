"""
    Implementation of the Vlasov-Poisson system using Gramian moment equations in 1D1D.
"""

"""
    Storage for the electric field values and related data for the Vlasov-Poisson system

    # Fields:
    - `E::Vector{Float64}`: Electric field values at spatial points
    - `Lx::Float64`: Length of the domain (x_max - x_min) for Poisson solve
    - `x_range::Vector{Float64}`: Spatial grid points
    - `initialized::Bool`: Flag to check if initialized
    - `ρ::Vector{Float64}`: Density values at spatial points (to be removed later)
    - `MP1::Int`: Number of moments + 1
    - `E_L2::Vector{Float64}`: Track L2 norm of electric field over time
    - `times::Vector{Float64}`: Track actual times when E_L2 is recorded
    - `variables`: Store solution variables at spatial points
    - `coordinates::Vector{Float64}`: Store spatial coordinates
"""
mutable struct ElectricFieldStorage
    E::Vector{Float64} # Electric field values at spatial points
    Lx::Float64 # Length of the domain (x_max - x_min) for Poisson solve
    x_range::Vector{Float64}    # Spatial grid points
    initialized::Bool   # Flag to check if initialized
    MP1::Int    # number of moments + 1
    E_L2::Vector{Float64}   # Track L2 norm of electric field over time
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
    variables#::Matrix{Float64} # Store solution variables at spatial points
    coordinates::Vector{Float64}    # Store spatial coordinates
end

# Global instance
const ELECTRIC_FIELD = ElectricFieldStorage(Float64[], 0.0, Float64[], false, 0, Float64[], Float64[],
    [Float64[]], Float64[]
)


"""
    Source term that applies the source term from the electric field to the moment equations
"""
function vlasov_poisson_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    E_field = linear_interpolation(ELECTRIC_FIELD.x_range, ELECTRIC_FIELD.E, extrapolation_bc = Interpolations.Line())
    
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
    return source
end

"""
    Solve the Poisson equation ∂E/∂x = ρ - ⟨ρ⟩ with periodic boundary conditions using FFT
"""
function solve_poisson_periodic_fft(ρ::AbstractVector{<:Real})
    n = length(ρ) # number of spatial points
    ρ̃  = ρ .- mean(ρ)  # neutralizing background
    ρk = fft(ρ̃)

    # build wavenumbers k consistent with FFT ordering
    # k = 0, 1, ..., floor(n/2), -ceil((n-1)/2), ..., -1
    k_int = [0:div(n,2); -div(n-1,2):-1]
    kx = (2π / ELECTRIC_FIELD.Lx) .* k_int

    # Solve for E in Fourier space
    Ek = similar(ρk)
    Ek[1] = 0 # k = 0 mode
    for j in 2:n
        if kx[j] != 0.0
            # Directly solves for E in Fourier space
            Ek[j] = ρk[j] / (im * kx[j])
        else
            Ek[j] = 0
        end
    end

    # Transform back to physical space
    E = real(ifft(Ek))
    return E
end

# Callback to solve Poisson equation globally at each timestep
"""
    Callback function to solve the Poisson equation and store electric field data
"""
function vlasov_poisson_callback(integrator)
    u = integrator.u
    t = integrator.t

    _, coordinates, variables = collect1dTreeArrays(integrator, cons2cons) # cons2cons only relevant for connectivity (first return argument) -> not relevant here
    if !ELECTRIC_FIELD.initialized
        ELECTRIC_FIELD.x_range = vec(coordinates[2:end-1])  # exclude ghost cells
        ELECTRIC_FIELD.initialized = true
    end

    # Extract density from solution variables
    # Helps for higher polynomial degrees, but could be optimized further
    ρ = variables[2:end-1, 1]  # First column is density 
    
    # Solve Poisson equation globally
    E = solve_poisson_periodic_fft(ρ)
    
    # Store the electric field
    ELECTRIC_FIELD.E = copy(E)

    # L2-norm of electric field
    # Energy = ||E(t,⋅)||_L2 = (∫ |E(t,x)|² dx)^(1/2)  (approximated via trapezoidal rule)
    E_L2 = trapz(ELECTRIC_FIELD.x_range, E.^2)^(1/2)
    push!(ELECTRIC_FIELD.E_L2, E_L2)
    push!(ELECTRIC_FIELD.times, t)  # Store the actual time
    
    return nothing
end

"""
    Create a DiscreteCallback for the Vlasov-Poisson system to solve Poisson equation at each timestep
"""
function vlasov_poisson_callback(;M, mesh, domain)
    # Reset the global storage to clear old data from previous runs
    empty!(ELECTRIC_FIELD.E)
    empty!(ELECTRIC_FIELD.E_L2)
    empty!(ELECTRIC_FIELD.times)
    ELECTRIC_FIELD.initialized = false
    ELECTRIC_FIELD.Lx = domain[end] - domain[1]
    
    # Set parameters for this run
    ELECTRIC_FIELD.MP1 = M+1
    
    return DiscreteCallback(
        (u, t, integrator) -> true,  # Always trigger at every step
        vlasov_poisson_callback,
        save_positions=(true, true), # ?(false, false),
        initialize = (c, u, t, integrator) -> vlasov_poisson_callback(integrator)  # Call at initialization
    )
end
