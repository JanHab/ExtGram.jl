mutable struct VlasovMaxwellStorage
    Ex::Vector{Float64} # Electric field values in x-direction
    Ey::Vector{Float64} # Electric field values in y-direction
    Ez::Vector{Float64} # Electric field values in z-direction
    Bx::Float64 # Magnetic field constant in x-direction
    By::Float64 # Magnetic field constant in y-direction
    Bz::Float64 # Magnetic field constant in z-direction
    Lx::Float64 # Length of the domain (x_max - x_min) for Poisson solve
    x_range::Vector{Float64}    # Spatial grid points
    initialized::Bool   # Flag to check if initialized
    MP1::Int    # number of moments + 1
    Ex_L2::Vector{Float64}   # Track L2 norm of electric field in x-direction over time
    Ey_L2::Vector{Float64}   # Track L2 norm of electric field in y-direction over time
    Ez_L2::Vector{Float64}   # Track L2 norm of electric field in z-direction over time
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
    variables#::Matrix{Float64} # Store solution variables at spatial points
    coordinates::Vector{Float64}    # Store spatial coordinates
end

# Global instance
const VLASOV_MAXWELL_FIELD = VlasovMaxwellStorage(Float64[], Float64[], Float64[], 0.0, 0.0, 0.0, 0.0, Float64[], false, 0, Float64[], Float64[], Float64[], Float64[], [Float64[]], Float64[]
)

"""
    moment_or_zero(u, equations, α, β, γ)

    `U_(αβγ)` if that moment is evolved by the system, and zero otherwise.

    The Vlasov source term references shifted multi-indices that can leave the system:
    negative powers (e.g. `β-1` for `β = 0`) and total order `M+1` (e.g. `α+1` for a moment
    of the order-M block). Neither contributes to the equations, so both are filtered to
    zero here rather than producing an empty index.
"""
@inline function moment_or_zero(u, equations::GramianMomentEquations1D3V, α::Int, β::Int, γ::Int)
    (α < 0 || β < 0 || γ < 0) && return zero(eltype(u))
    (α + β + γ > equations.M) && return zero(eltype(u))
    p = @inbounds equations._pos[α + 1, β + 1, γ + 1]
    return p == 0 ? zero(eltype(u)) : @inbounds u[p]
end

function vlasov_maxwell_source(u, x, t, equations::GramianMomentEquations1D3V{Mp1, N, MC, NC, NS}) where {Mp1, N, MC, NC, NS}
    # ToDo: Move this out of the source term to only call it once, not multiple times
    Ex_field = linear_interpolation(VLASOV_MAXWELL_FIELD.x_range, VLASOV_MAXWELL_FIELD.Ex, extrapolation_bc = Interpolations.Line())
    # Ignore for the moment, as it is set to zero
    # Ey_field = linear_interpolation(VLASOV_MAXWELL_FIELD.x_range, VLASOV_MAXWELL_FIELD.Ey, extrapolation_bc = Interpolations.Line())
    # Ez_field = linear_interpolation(VLASOV_MAXWELL_FIELD.x_range, VLASOV_MAXWELL_FIELD.Ez, extrapolation_bc = Interpolations.Line())

    Ex_local = Ex_field(x[1])
    # Ey/Ez are set to zero for the moment (see the commented interpolants above)
    Ey_local = zero(Ex_local)
    Ez_local = zero(Ex_local)

    source = MVector{Mp1, Float64}(undef)
    source[1] = 0.0                    # No source for mass conservation

    # higher moments
    @inbounds for k in 2:Mp1
        # Powers (l, m, n) of the moment this equation evolves
        α, β, γ = equations._pow[k]

        # l * E_x * u_{l-1, m, n} + m * E_y * u_{l, m-1, n} + n * E_z * u_{l, m, n-1}
        source[k] = α * Ex_local * moment_or_zero(u, equations, α - 1, β, γ) +
                    β * Ey_local * moment_or_zero(u, equations, α, β - 1, γ) +
                    γ * Ez_local * moment_or_zero(u, equations, α, β, γ - 1)

        # Magnetic field contributions, from (v x B) . grad_v:
        #   B_x: m u_{l, m-1, n+1} - n u_{l, m+1, n-1}
        #   B_y: n u_{l+1, m, n-1} - l u_{l-1, m, n+1}
        #   B_z: l u_{l-1, m+1, n} - m u_{l+1, m-1, n}
        source[k] += (
            α * (
                VLASOV_MAXWELL_FIELD.Bz * moment_or_zero(u, equations, α - 1, β + 1, γ) - VLASOV_MAXWELL_FIELD.By * moment_or_zero(u, equations, α - 1, β, γ + 1)
            ) + β * (
                VLASOV_MAXWELL_FIELD.Bx * moment_or_zero(u, equations, α, β - 1, γ + 1) - VLASOV_MAXWELL_FIELD.Bz * moment_or_zero(u, equations, α + 1, β - 1, γ)
            ) + γ * (
                VLASOV_MAXWELL_FIELD.By * moment_or_zero(u, equations, α + 1, β, γ - 1) - VLASOV_MAXWELL_FIELD.Bx * moment_or_zero(u, equations, α, β + 1, γ - 1)
            )
        )
    end
    
    return source
end

# Callback to solve Poisson equation globally at each timestep
"""
    Callback function to solve the Vlasov-Maxwell equation and store electric field data
"""
function update_electric_field!(semi, u_ode)
    _, coordinates, variables = collect1dTreeArrays(semi, u_ode, cons2cons)
    if !VLASOV_MAXWELL_FIELD.initialized
        VLASOV_MAXWELL_FIELD.x_range = vec(coordinates[1:end-1])
        VLASOV_MAXWELL_FIELD.initialized = true
    end
    VLASOV_MAXWELL_FIELD.Ex = solve_poisson_periodic_fft(variables[1:end-1, 1], VLASOV_MAXWELL_FIELD)
    return nothing
end

function rhs_vlasov_maxwell!(du_ode, u_ode, semi, t)
    update_electric_field!(semi, u_ode)   # E belongs to *this* stage's state
    Trixi.rhs!(du_ode, u_ode, semi, t)
    return nothing
end

function vlasov_maxwell_diagnostics(integrator)
    update_electric_field!(integrator.p, integrator.u)  # E for u(t_n), not the last stage
    Ex = VLASOV_MAXWELL_FIELD.Ex
    push!(VLASOV_MAXWELL_FIELD.Ex_L2, trapz(VLASOV_MAXWELL_FIELD.x_range, Ex.^2)^(1/2))
    push!(VLASOV_MAXWELL_FIELD.times, integrator.t)
    return nothing
end

"""
    Create a DiscreteCallback for the Vlasov-Maxwell system to solve the Vlasov-Maxwell equation at each timestep
"""
function vlasov_maxwell_callback(;Bx, By, Bz, M, mesh, domain)
    # Reset the global storage to clear old data from previous runs
    empty!(VLASOV_MAXWELL_FIELD.Ex)
    empty!(VLASOV_MAXWELL_FIELD.Ey)
    empty!(VLASOV_MAXWELL_FIELD.Ez)
    empty!(VLASOV_MAXWELL_FIELD.Ex_L2)
    empty!(VLASOV_MAXWELL_FIELD.Ey_L2)
    empty!(VLASOV_MAXWELL_FIELD.Ez_L2)
    empty!(VLASOV_MAXWELL_FIELD.times)
    VLASOV_MAXWELL_FIELD.initialized = false
    VLASOV_MAXWELL_FIELD.Lx = domain[end] - domain[1]
    VLASOV_MAXWELL_FIELD.Bx = Bx
    VLASOV_MAXWELL_FIELD.By = By
    VLASOV_MAXWELL_FIELD.Bz = Bz
    
    # Set parameters for this run
    VLASOV_MAXWELL_FIELD.MP1 = M+1
    
    return DiscreteCallback(
        (u, t, integrator) -> true,  # Always trigger at every step
        vlasov_maxwell_diagnostics,
        save_positions=(false, false),  # the field solve lives in rhs_vlasov_maxwell!, nothing to save around
        initialize = (c, u, t, integrator) -> vlasov_maxwell_diagnostics(integrator)  # seeds E(u_0) and records t = 0
    )
end

struct InitialConditionsVlasovMaxwellLandauDamping{N}
    v1::Float64
    v2::Float64
    v3::Float64
    α::Float64
    k::Float64
    # Moments of the Maxwellian for ρ = 1, in the ordering of the evolved moments.
    # Only ρ varies with x here (v and θ are uniform) and the moments are linear in ρ,
    # so the state at a node is just this vector scaled by ρ(x).
    unit_moments::SVector{N, Float64}
    equations::GramianMomentEquations1D3V

    function InitialConditionsVlasovMaxwellLandauDamping(v1::Float64, v2::Float64, v3::Float64, α::Float64, k::Float64, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
        M = equations.M
        full = multi_index_list(M)
        position_in_full = Dict(ix => i for (i, ix) in enumerate(full))
        gather = [position_in_full[ix] for ix in equations._U_t_index]

        unit_moments = convective_moments_1D3D(M, 1.0, (v1, v2, v3), 1.0)[gather]
        @assert length(unit_moments) == Mp1 == equations.N_equations

        return new{Mp1}(v1, v2, v3, α, k, SVector{Mp1, Float64}(unit_moments), equations)
    end
end

function (ic::InitialConditionsVlasovMaxwellLandauDamping)(coords, t, equations::GramianMomentEquations1D3V)
    ρx = 1 + ic.α * cos(ic.k * coords[1])
    return ρx * ic.unit_moments
end
