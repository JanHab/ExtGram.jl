struct BGKEquations1D{N} <: Trixi.AbstractEquations{1, N}
    N::Int  # Number of equations to be solved
    c_l::Real   # Lower bound of the velocity domain
    c_u::Real   # Upper bound of the velocity domain
    c_vec::SVector{N, Float64}  # Velocity grid points
    Kn::Real    # Knudsen number
    function BGKEquations1D(N::Integer, c_l::Real, c_u::Real, Kn::Real)
        @assert N > 1 "Number of equations N must be positive"
        @assert c_l < c_u "Lower bound c_l must be smaller than upper bound c_u"
        c_vec = SVector{N, Float64}(Tuple(LinRange(c_l, c_u, N)))
        new{N}(N, c_l, c_u, c_vec, Kn)
    end
end
dc(equations::BGKEquations1D) = (equations.c_u - equations.c_l) / (equations.N - 1)

# * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), ::BGKEquations1D{N}) where {N} = ntuple(i->"f^{$(i-1)}", N)

function Trixi.flux(f, orientation::Integer, equations::BGKEquations1D{N}) where {N}
    return SVector(
        ntuple(i -> f[i] * equations.c_vec[i], N)
    )
end

# todo:
# Overloading for only one relaxation_source function
function relaxation_source(f, x, t, equations::BGKEquations1D{N}) where {N} 
    ρ = trapz(equations.c_vec, f) # density
    v = trapz(equations.c_vec, f .* equations.c_vec) / ρ # velocity
    θ = trapz(equations.c_vec, f .* (equations.c_vec .- v).^ 2) / ρ # temperature

    f_Maxwellian = Maxwellian(ρ, v, θ).(equations.c_vec)
    # @assert all(f_Maxwellian .>= 0.0)
    return (-1.0 / equations.Kn) * (f .- f_Maxwellian)

    # Assemble weight matrix A
    A = zeros(3, N)
    A[1, :] .= dc(equations) # ρ = Σ w_i f_i -> A[1, :] = w_i = dc
    A[2, :] .= dc(equations) .* equations.c_vec # ρ v = Σ w_i c_i f_i -> A[2, :] = w_i * c_i
    A[3, :] .= dc(equations) .* equations.c_vec.^2 # ρ e = Σ w_i c_i^2 f_i -> A[3, :] = w_i * c_i^2

    # Assemble left-hand side of Lagrange-multiplier system
    Ã = zeros(N+3, N+3)
    Ã[1:N, 1:N] .= I(N) # identity matrix
    Ã[N+1:N+3, 1:N] .= A
    Ã[1:N, N+1:N+3] .= A'
    # Ã[N+1:N+3, N+1:N+3] .= # zeros

    # Assemble right hand side vector b
    b = zeros(N+3)
    # b[1:N] = 0
    r = SVector{3}(ρ, ρ .* v, ρ .* θ)
    Δr = r - A * f_Maxwellian
    b[N+1:N+3] .= Δr - A * f_Maxwellian

    # Solve for Δf
    # println("Solving BGK Lagrange-multiplier system...")
    Δf_full = Ã \ b
    # println("...done.")
    Δf = similar(f)
    Δf .= Δf_full[1:N] # updates for distribution function
    # println("Δf = ", Δf)
    λ = Δf_full[N+1:N+3] # Lagrange multipliers
    # println("λ = ", λ)

    Δf = A' * inv(A * A') * Δr

    f̃ = f_Maxwellian + Δf
    # println("f̃ = ", f̃)

    ρ̃ = trapz(equations.c_vec, f̃)
    # println("ρ: ", ρ)
    # println("ρ̃: ", ρ̃)
    # println("abs(ρ - ρ̃) / ρ = ", abs(ρ - ρ̃) / ρ)

    return (-1.0 / equations.Kn) * (f .- f̃)
end
zero_source(f, x, t, equations::BGKEquations1D{N}) where {N} = SVector{N}(ntuple(i->0.0, N))

function flux_jacobian(u, equations::BGKEquations1D{N}) where {N}
    return Diagonal(equations.c_vec)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
function Trixi.max_abs_speed_naive(f_l, f_r, orientation::Integer, equations::BGKEquations1D)
    λ_l = Trixi.max_abs_speeds(f_l, equations)
    λ_r = Trixi.max_abs_speeds(f_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end

function Trixi.max_abs_speeds(f, equations::BGKEquations1D)
    # estimate the flux Jacobian eigenvalues by means of Gerschgorin
    # return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
    return maximum(abs.(real.(eigen(flux_jacobian(f, equations)).values)))
end


# ========================================================================================== #

# Initial Conditions
struct InitialConditionsBGK{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsBGK(f_left, f_right, equations::BGKEquations1D{N}) where {N}
        left = f_left.(equations.c_vec) #? ones(SVector{N}) * f_left
        right = f_right.(equations.c_vec) #? ones(SVector{N}) * f_right
        @assert all(left .>= 0.0) && all(right .>= 0.0)

        return new{N}(left, right)
    end
end

function (ic::InitialConditionsBGK)(coords, t, equations::BGKEquations1D{N}) where {N}
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end




# ========================================================================================== #

# physical variables
# u_i = ∫ c^i f(c) dc, i=1,...,N
# u_i(x) = ∫ c_k^i f(c_k)|_x dc, k=1,...,N
function ρ_v_θ_p_BGK(f::Vector, equations::BGKEquations1D)
    L = length(f)
    Nloc = L ÷ equations.N
    Fmat = reshape(f, equations.N, Nloc)
    ρ = dc(equations) * sum(Fmat, dims=1)
    v = dc(equations) * sum(Fmat .* equations.c_vec, dims=1) ./ ρ
    Θ = dc(equations) * sum(Fmat .* (equations.c_vec .- v).^ 2, dims=1) ./ ρ
    p = ρ .* Θ
    return vec(ρ), vec(v), vec(Θ), vec(p)
end

function plot_ρ_v_p_bgk(ρ, v, p, x; xlims=(-2.0, 2.0))
    # Plot physical variables
    # create primary plot and plot density and pressure on it
    plt = plot(x, ρ;
        xlabel = "x",
        label = "ρ",
        color = :blue,
        legend = :topleft,
    )

    # plot pressure on the same (primary) axis
    plot!(plt, x, p;
        label = "p",
        color = :green,
    )

    # create a twin y-axis for velocity
    ax2 = twinx(plt)
    plot!(ax2, x, v;
        label = "v",
        color = :red,
        legend = :topright,
    )

    # set x limits on the primary axis (affects both)
    xlims!(plt, xlims)

    return plt
end


# setup for 1D Riemann
function setupBGK1DRiemann(
    N, Kn,
    c_l, c_u, 
    ic_left, ic_right;
    base_tree_level = 8,
    surface_flux = flux_lax_friedrichs,
    volume_flux = flux_central,
    polydeg = 1,        # DG polynomial degree
    domain = (-5.0, 5.0),
    )
    equations = BGKEquations1D(N, c_l, c_u, Kn)
    initial_condition = InitialConditionsBGK(
        ic_left,
        ic_right,
        equations
    )

    #= set up semidiscretization =#
    basis = LobattoLegendreBasis(polydeg)

    # ==================================== volume integral treatment ==================================== #

    # shock capturing
    # indicator_sc = IndicatorHennemannGassner(
    #     equations, 
    #     basis,
    #     alpha_max = 1.0,
    #     alpha_min = 0.01,
    #     alpha_smooth = true,
    #     variable = (u, eqns)->u[1]*u[3]
    # ) # ? Which variable?

    # volume_integral = VolumeIntegralShockCapturingHG(
    #     indicator_sc;
    #     volume_flux_dg = volume_flux,
    #     volume_flux_fv = surface_flux
    # )
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

    # =================================================================================================== #


    solver = DGSEM(basis, surface_flux, volume_integral)

    mesh = TreeMesh(
        (domain[1],), 
        (domain[2],), 
        initial_refinement_level=base_tree_level, 
        n_cells_max=10_000, 
        periodicity=false
    )

    boundary_conditions = (
        x_neg = BoundaryConditionDirichlet(initial_condition), 
        x_pos = BoundaryConditionDirichlet(initial_condition)
    )


    return basis, mesh, equations, initial_condition, solver, boundary_conditions
end




# ========================================================================================== #
# Vlasov-Poisson specific callbacks and initial conditions

# Global storage for the electric field
# This will be updated during each RHS evaluation
mutable struct ElectricFieldStorageBGK
    E::Vector{Float64}
    domain::Tuple{Float64, Float64}
    x_range::Vector{Float64}
    n::Int
    initialized::Bool
    ρ::Vector{Float64}
    N::Int
    E_L2::Vector{Float64}
    times::Vector{Float64}  # Track actual times when E_L2 is recorded
    variables#::Matrix{Float64}
    c_vec::Vector{Float64}  # Velocity grid points
    dc::Float64
end

# Global instance
const ELECTRIC_FIELD_BGK = ElectricFieldStorageBGK(Float64[], (0.0, 0.0), Float64[], 0, false, Float64[], 0, Float64[], Float64[],
    [Float64[]], # todo: remove this later
    Float64[],
    0.0
)

# Source term that solves Poisson globally and applies local source
function vlasov_poisson_source(f, x, t, equations::BGKEquations1D{N}) where {N}
    # todo: hard-coded for the moment
    x_range = ELECTRIC_FIELD_BGK.x_range
    E_field = linear_interpolation(x_range, ELECTRIC_FIELD_BGK.E, extrapolation_bc = Interpolations.Line())
    
    # Evaluate electric field at position x
    E_local = E_field(x[1])

    # Apply source terms
    source = MVector{N, Float64}(undef)

    ∂ᵥf = similar(f)
    dc = ELECTRIC_FIELD_BGK.dc
    # central differences in velocity space
    for i in 1:N
        if i == 1
            ∂ᵥf[i] = (f[i+1] - f[i]) / dc
        elseif i == N
            ∂ᵥf[i] = (f[i] - f[i-1]) / dc
        else
            ∂ᵥf[i] = (f[i+1] - f[i-1]) / (2.0 * dc)
        end
    end
    source = -E_local .* ∂ᵥf

    # return positive value of source, as it is on the lhs but with a minus -> positive source term
    return source
end

function solve_poisson_periodic_fft_BGK(ρ::AbstractVector{<:Real})
    n = length(ρ)
    Lx = ELECTRIC_FIELD_BGK.domain[end] - ELECTRIC_FIELD_BGK.domain[1]
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
function vlasov_poisson_callback_BGK(integrator)
    u = integrator.u
    t = integrator.t

    connectivity, coordinates, variables = collect1dTreeArrays(integrator, cons2cons) # cons2cons only relevant for connectivity -> not relevant here

    Fmat = variables[2:end-1, :]  # First column is density # todo: why 2:end-1? What's wrong here? The first and last entry seem to be off, though.
    ρ = vec(ELECTRIC_FIELD_BGK.dc * sum(Fmat, dims=2))
    ELECTRIC_FIELD_BGK.variables = variables # todo: remove this later
    ELECTRIC_FIELD_BGK.ρ = ρ
    ELECTRIC_FIELD_BGK.x_range = vec(coordinates[2:end-1])

    # Solve Poisson equation globally
    E = solve_poisson_periodic_fft_BGK(ρ)
    # Store the electric field
    ELECTRIC_FIELD_BGK.E = copy(E)

    # Add energy vector
    # Energy = ||E(t,⋅)||_L2 = (∫ |E(t,x)|² dx)^(1/2)  (approximated via trapezoidal rule) -> L2-norm of electric field
    E_L2 = trapz(ELECTRIC_FIELD_BGK.x_range, E.^2)^(1/2)
    push!(ELECTRIC_FIELD_BGK.E_L2, E_L2)
    push!(ELECTRIC_FIELD_BGK.times, t)  # Store the actual time

    ELECTRIC_FIELD_BGK.initialized = true
    
    return nothing
end

# Create the callback - triggers after each iteration
function vlasov_poisson_callback_BGK(;N, mesh, domain, equations)
    # Reset the global storage to clear old data from previous runs
    empty!(ELECTRIC_FIELD_BGK.E)
    empty!(ELECTRIC_FIELD_BGK.E_L2)
    empty!(ELECTRIC_FIELD_BGK.times)
    empty!(ELECTRIC_FIELD_BGK.ρ)
    ELECTRIC_FIELD_BGK.n = 0
    ELECTRIC_FIELD_BGK.initialized = false
    ELECTRIC_FIELD_BGK.domain = domain
    ELECTRIC_FIELD_BGK.c_vec = equations.c_vec
    ELECTRIC_FIELD_BGK.dc = dc(equations)
    
    # Set parameters for this run
    ELECTRIC_FIELD_BGK.N = N
    
    return DiscreteCallback(
        (u, t, integrator) -> true,  # Always trigger at every step
        vlasov_poisson_callback_BGK,
        save_positions=(true, true), # ?(false, false),
        initialize = (c, u, t, integrator) -> vlasov_poisson_callback_BGK(integrator)  # Call at initialization
    )
end

struct InitialConditionsCosine_BGK{N}
    ρ0::Float64
    ϵ::Float64
    v0::Float64
    θ0::Float64
    k::Float64

    function InitialConditionsCosine_BGK(ρ0::Float64, ϵ::Float64, v0::Float64, θ0::Float64, k::Float64, eqns::BGKEquations1D{N}) where {N}
        return new{N}(ρ0, ϵ, v0, θ0, k)
    end
end

function (ic::InitialConditionsCosine_BGK)(coords, t, equations::BGKEquations1D{N}) where {N}
    ρx = ic.ρ0 * (1 + ic.ϵ * cos(ic.k * coords[1]))
    f = Maxwellian(ρx, ic.v0, ic.θ0)
    # return convective_moments(f, Val(Mp1))
    return f.(equations.c_vec)
end
