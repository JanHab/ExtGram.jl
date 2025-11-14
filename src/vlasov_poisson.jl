# Global storage for the electric field
# This will be updated during each RHS evaluation
mutable struct ElectricFieldStorage
    E::Vector{Float64} # Electric field values at spatial points
    Lx::Float64 # Length of the domain (x_max - x_min) for Poisson solve
    x_range::Vector{Float64}    # Spatial grid points
    initialized::Bool   # Flag to check if initialized
    ρ::Vector{Float64} # todo: remove this later
    MP1::Int    # number of moments + 1
    E_L2::Vector{Float64}   # Track L2 norm of electric field over time
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
    variables#::Matrix{Float64} # Store solution variables at spatial points
    coordinates::Vector{Float64}    # Store spatial coordinates
end

# Global instance
const ELECTRIC_FIELD = ElectricFieldStorage(Float64[], 0.0, Float64[], false, Float64[], 0, Float64[], Float64[],
    [Float64[]], Float64[]
)

# Source term that solves Poisson globally and applies local source
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

function solve_poisson_periodic_fft(ρ::AbstractVector{<:Real})
    n = length(ρ) # number of spatial points
    ρ̃  = ρ .- mean(ρ) #!ρ .- 1#! mean(ρ)  # neutralizing background
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
    ρ = variables[2:end-1, 1]  # First column is density # todo: why 2:end-1? What's wrong here? The first and last entry seem to be off, though.

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

# Create the callback - triggers after each iteration
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

struct InitialConditionsLandauDamping{N}
    ρ0::Float64
    ϵ::Float64
    v0::Float64
    θ0::Float64
    k::Float64

    function InitialConditionsLandauDamping(ρ0::Float64, ϵ::Float64, v0::Float64, θ0::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(ρ0, ϵ, v0, θ0, k)
    end
end

function (ic::InitialConditionsLandauDamping)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    ρx = ic.ρ0 * (1 + ic.ϵ * cos(ic.k * coords[1]))
    f = Maxwellian(ρx, ic.v0, ic.θ0)
    return convective_moments(f, Val(Mp1))
end

struct InitialConditionsTwoStream{N}
    ϵ::Float64
    k::Float64

    function InitialConditionsTwoStream(ϵ::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(ϵ, k)
    end
end

# function (ic::InitialConditionsTwoStream)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
#     v0 = 2.0
#     f = c -> 0.5/sqrt(2π) * (exp.(-(c - v0).^2 ./ 2) .+ exp.(-(c + v0).^2 ./ 2)) .* (1 .+ ic.ϵ * cos(ic.k * coords[1]))
#     ξ, w = gausshermite(Mp1+1)
#     C = sqrt(2.0) .* ξ
#     fw = f.(C) .* w .* exp.(ξ.^2) * sqrt(2.0)
#     return SVector{Mp1,Float64}(ntuple(n->sum(C .^(n-1) .* fw), Mp1))
# end


function (ic::InitialConditionsTwoStream)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    max(c) = 1/sqrt(2*π) * exp(-c^2 / 2) * c^2 .* (1 + ic.ϵ * cos(ic.k * coords[1])) # note: normalized in 1D velocity space 
    # copying from convective_moments for standard Maxwellian
    ξ, w = gausshermite(Mp1+1) # +1 for good measure, should not be necessary
    C = sqrt(2*1.0) .* ξ; c = C .+ 0.0
    fw = max.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*1.0)
    return SVector{Mp1,Float64}(ntuple(n->sum(c .^(n-1) .* fw), Mp1))
end

# using QuadGK

# struct InitialConditionsTwoStream{N}
#     ϵ::Float64
#     k::Float64
#     v0::Float64

#     # default constructor keeps the old call signature InitialConditionsTwoStream(ϵ,k,eqns)
#     function InitialConditionsTwoStream(ϵ::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
#         return new{Mp1}(ϵ, k, 2.5) # default beam speed v0=2.5 (in thermal units)
#     end

#     # alternative constructor allowing explicit beam speed
#     function InitialConditionsTwoStream(ϵ::Float64, k::Float64, v0::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
#         return new{Mp1}(ϵ, k, v0)
#     end
# end

# function (ic::InitialConditionsTwoStream)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
#     x = coords[1]
#     # density modulation in x
#     ρscale = 1.0 * (1.0 + ic.ϵ * cos(ic.k * x))

#     # two symmetric beams shifted by ±v0 (each with half the density)
#     v0 = ic.v0
#     # build the total distribution f(c) = ρscale * 0.5*(M(v0) + M(-v0)) with thermal θ=1.0
#     Mplus  = Maxwellian(1.0, v0, 1.0)
#     Mminus = Maxwellian(1.0, -v0, 1.0)

#     f_total(c) = ρscale * 0.5 * (Mplus(c) + Mminus(c))

#     # compute moments by numerical integration over velocity (use quadgk on (-Inf, Inf))
#     moments = ntuple(n -> quadgk(c -> c^(n-1) * f_total(c), -Inf, Inf)[1], Mp1)
#     return SVector{Mp1,Float64}(moments)
# end