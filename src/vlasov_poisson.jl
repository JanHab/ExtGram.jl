# Global storage for the electric field
# This will be updated during each RHS evaluation
mutable struct ElectricFieldStorage
    E::Vector{Float64}
    domain::Tuple{Float64, Float64}
    n::Int
    initialized::Bool
    ρ::Vector{Float64} # todo: remove this later
    MP1::Int
    E_L2::Vector{Float64}
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
end

# Global instance
const ELECTRIC_FIELD = ElectricFieldStorage(Float64[], (-0.0, 0.0), 0, false, Float64[], 0, Float64[], Float64[])

function solve_poisson_periodic_fft(ρ::AbstractVector{<:Real}, domain::Tuple{Float64,Float64})
    n = length(ρ)
    Lx = domain[2] - domain[1]
    ρ̃ = ρ .- 1#! mean(ρ)  # neutralizing background
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

# Source term that solves Poisson globally and applies local source
function vlasov_poisson_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    # todo: hard-coded for the moment
    x_min, x_max = ELECTRIC_FIELD.domain
    x_range = range(x_min, x_max, length=ELECTRIC_FIELD.n)
    E_field = interpolate((x_range,), ELECTRIC_FIELD.E, Gridded(Linear()))
    
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

# Callback to solve Poisson equation globally at each timestep
function vlasov_poisson_callback(integrator)
    u = integrator.u
    t = integrator.t
    
    # Extract density from the full solution vector
    MP1 = ELECTRIC_FIELD.MP1
    n_cells = length(u) ÷ MP1
    ELECTRIC_FIELD.n = n_cells
    
    # Reshape to extract density
    U_matrix = reshape(u, MP1, n_cells)
    ρ = U_matrix[1, :]  # First row is density
    ELECTRIC_FIELD.ρ = ρ

    # Solve Poisson equation globally
    E = solve_poisson_periodic_fft(ρ, ELECTRIC_FIELD.domain)
    # Store the electric field
    ELECTRIC_FIELD.E = copy(E)

    # Add energy vector
    # Energy = ||E(t,⋅)||_L2 = (∫ |E(t,x)|² dx)^(1/2)  (approximated via trapezoidal rule)
    # L2-norm of electric field
    Δx = (ELECTRIC_FIELD.domain[2] - ELECTRIC_FIELD.domain[1]) / n_cells
    E_L2 = (sum(E.^2) * Δx)^(1/2)
    # E_L2 = 1/2 * sum(E.^2) * (ELECTRIC_FIELD.domain[2] - ELECTRIC_FIELD.domain[1]) / n_cells
    push!(ELECTRIC_FIELD.E_L2, E_L2)
    push!(ELECTRIC_FIELD.times, t)  # Store the actual time

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
