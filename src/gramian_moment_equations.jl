
# (note: replacing all zeroes() calls with corresponding SVector and MVector calls would speed up things at the cost of arbitrary vector input types)

# * Basic struct for equations
"""
struct GramianMomentEquations1D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}

A struct representing one-dimensional Gramian moment equations for kinetic theory.

    closure::Symbol
- `inv_Kn::RealT`: The inverse of the Knudsen number, used in the relaxation source term
    function GramianMomentEquations1D(M::Integer, Knudsen::Real, closure::Union{Symbol,String}=:ExtGram; χ_set="optimal")
- `n::Int`: Internal parameter derived from the number of moments
- `extended::Bool`: Flag indicating whether to use extended formulation

# Constructor
    GramianMomentEquations1D(M::Integer, Knudsen::Real, extended=true::Bool)

Constructs a `GramianMomentEquations1D` system with `M` moments and specified Knudsen number.

# Arguments
- `M::Integer`: Number of moments (must be greater than 1)
- `Knudsen::Real`: Knudsen number for the system
- `closure::Union{Symbol,String}`: Closure type, can be "Gram", "ExtGram", or "Grad"
# Notes
- The system has `M+1` variables in total
- The parameter `χ` is calculated differently depending on whether `M` is even or odd
"""
struct GramianMomentEquations1D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}
    inv_Kn::RealT # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    n::Int
    closure::Symbol

    function GramianMomentEquations1D(M::Integer, Knudsen::Real, closure::String="ExtGram"; χ_set="optimal")
        @assert M > 1
        if iseven(M)
            n = Int(M/2)
            χ = (n+1)/n
        else
            n = Int((M+1)/2)
            χ = (n+1)/(2n)
        end
        # overwrite if χ_set is given explicitly, only relevant for extended Gramian closure
        if χ_set !== "optimal"
            χ = χ_set
        end
        closure_value = :ExtGram # default option
        if closure == "Gram"
            closure_value = iseven(M) ? :GramEven : :GramOdd
        elseif closure == "ExtGram"
            closure_value = iseven(M) ? :ExtGramEven : :ExtGramOdd
        elseif closure == "Grad"
            closure_value = :Grad
        elseif closure == "MaxEnt"
            closure_value = :MaxEnt
        else
            error("Unknown closure type: $closure. Supported types are \"Gram\", \"ExtGram\", \"Grad\" and \"MaxEnt (work in progress)\n.")
        end
        new{M+1, n, typeof(Knudsen)}(inv(Knudsen), χ, n, closure_value)
    end
end

# * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), ::GramianMomentEquations1D{Mp1}) where {Mp1} = ntuple(i->"u^($(i-1))", Mp1)
Trixi.varnames(::typeof(cons2prim), ::GramianMomentEquations1D{Mp1}) where {Mp1} = ntuple(i->"w^($(i-1))", Mp1)

# * basic variable definitions to be used as e.g. sc indicator variable
Trixi.density(u, eqns::GramianMomentEquations1D{Mp1}) where {Mp1} = u[1]



"""
    Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D{Mp1}) where {Mp1}

    Flux function: F(U) = u_{k+1}, k=0,…,M where u_{M+1} = C(u0, u1, …, uM)
"""
function Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    # First MP1-1 flux components from shifed moments
    # ∂_t u^k + \partial_x u^{k+1} = ... (0)
    # Last u from closure
    # return SVector(ntuple(i->u[i+1], Mp1-1)..., closure(u, equations))
    return SVector(ntuple(i->u[i+1], Mp1-1)..., closure(u, equations, Val(equations.closure)))
end
isinvertible(A::Matrix{Float64}) = !isapprox(det(BigFloat.(A)), 0, atol = 1e-18)


########################## Closure ##########################

# Overloading to determine Gram/ExtGram/Grad and even/odd case
"""
    closure(u, equations::GramianMomentEquations1D)

    Decides which closure implementation is used based on the equations.closure value
"""
closure(u, equations::GramianMomentEquations1D) = closure(u, equations, Val(equations.closure))

# -------------------------
# Gramian closure 
# -------------------------
# Decides which closure implementation is used based on the number of Moments (even vs. odd)
# closure(u, equations::GramianMomentEquations1D, ::Val{:Gram}) = closure(u, equations, Val(iseven(length(u)-1) ? :GramEven : :GramOdd))
# # -------------------------
# closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGram}) = closure(u, equations, Val(iseven(length(u)-1) ? :ExtGramEven : :ExtGramOdd))

# Even case (classical and extended)
"""
    closure(u, equations::GramianMomentEquations1D, ::Val{true}; verbose_=false)

    Closure for the even case
"""
function closure(u, equations::GramianMomentEquations1D, ::Val{:GramEven})
    M = length(u)-1 # u[0, ..., M]
    @assert iseven(M)
    n = equations.n 
    invG_nm1 = inv(gramian(u, n-1))
    return @views(
        u[n+2:2n+1]'*invG_nm1*u[n+1:2n]
    )
end

function closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGramEven})
    M = length(u)-1 # u[0, ..., M]
    @assert iseven(M)
    n = equations.n 
    invG_nm1 = inv(gramian(u, n-1))
    invG_nm2 = inv(gramian(u, n-2))
    return @views (
        u[n+2:2n+1]'*invG_nm1*u[n+1:2n]
    ) + equations.χ * (
        u[2n+1] - u[n+1:2n]'*invG_nm1*u[n+1:2n]
    ) / (
        u[2n-1] - u[n:2n-2]'*invG_nm2*u[n:2n-2]
    ) * (
        u[2n] - u[n+1:2n-1]'*invG_nm2*u[n:2n-2]
    )
end

# Odd case (classical and extended)
"""
    closure(u, equations::GramianMomentEquations1D, ::Val{false})

    Closure for the odd case
"""
function closure(u, equations::GramianMomentEquations1D, ::Val{:GramOdd})
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    invG_nm1 = inv(gramian(u, n-1))
    return @views (
        u[n+1:2n]'*invG_nm1*u[n+1:2n]
    )
end

function closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGramOdd})
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    if n==2 # M=3
        invG_1 = inv(gramian(u, 1))
        invG_0 = 1 / u[1]
        # invG_2 = inv(gramian(u, 2))
        # Index mapping from SymPy (0-based): u0→u[1], u1→u[2], u2→u[3], u3→u[4], u4→u[5], u5→u[6]
        closure = @views (
            u[3:4]'*invG_1*u[3:4]
        ) + (
            u[3] - u[2]*invG_0*u[2]
        ) * (
            2 + (
                u[4] - u[3]*invG_0*u[2]
            )^2 / (
                u[3] - u[2]*invG_0*u[2]
            )^2
        )
    elseif n==3 # M=5
        invG_1 = inv(gramian(u, 1))
        invG_0 = 1 / u[1]
        invG_2 = inv(gramian(u, 2))
        A = @views (
            1 + (
                u[5] - u[3:4]'*invG_1*u[3:4]
            ) / (
                u[3] - u[2] * invG_0 * u[2]
            ) + 0.5 * (
                (
                    (
                        u[6] - u[4:5]'*invG_1*u[3:4]
                    ) / (
                        u[5] - u[3:4]'*invG_1*u[3:4]
                    ) - (
                        u[4] - u[3] * invG_0 * u[2]
                    ) / (
                        u[3] - u[2] * invG_0 * u[2]
                    )
                )^2 + (
                    (
                        u[6] - u[4:5]'*invG_1*u[3:4]
                    ) / (
                        u[5] - u[3:4]'*invG_1*u[3:4]
                    ) - (
                        2 * (
                            u[4] - u[3] * invG_0 * u[2]
                        ) / (
                            u[3] - u[2] * invG_0 * u[2]
                        )
                    )
                )^2
            )
        )
        
        closure = @views (
            u[4:6]'*invG_2*u[4:6]
        ) + (
            u[5] - u[3:4]'*invG_1*u[3:4]
        ) * A
    else
        sigma_nn(u, n) = u[2n+1] - u[n+1:2n]'*inv(gramian(u, n-1))*u[n+1:2n]
        sigma_nmnm(u, n) = u[2n-1] - u[n:2n-2]'*inv(gramian(u, n-2))*u[n:2n-2]
        sigma_nnp(u, n) = u[2n+2] - u[n+2:2n+1]'*inv(gramian(u, n-1))*u[n+1:2n]
        sigma_nmn(u, n) = u[2n] - u[n+1:2n-1]'*inv(gramian(u, n-2))*u[n:2n-2]

        an(u, n) = sigma_nnp(u, n) / sigma_nn(u, n) - sigma_nmn(u, n) / sigma_nmnm(u, n)
        bn(u, n) = sigma_nn(u, n) / sigma_nmnm(u, n)


        bk_sum = 1.0
        anm = an(u, n-1)
        ak_sum = (
            anm
        )^2 + (
            anm - (
                u[4] - u[3] / u[1] * u[2]
            ) / (
                u[3] - u[2] / u[1] * u[2]
            )
        )^2# (anm - a0)^2 + (anm - a1)^2
        for k=2:n-1
            bk_sum += bn(u, k)
            ak_sum += (anm - an(u, k))^2
        end
        closure = @views(
            u[n+1:2n]' * inv(gramian(u, n-1)) * u[n+1:2n]
        ) + (
            u[2n-1] - u[n:2n-2]'*inv(gramian(u, n-2))*u[n:2n-2]
        ) * (
            4 / (2*(n-1)) * bk_sum + 2 / (2*(n-1)) * ak_sum
        )
    end

    return closure
end


# -------------------------
# Grad closure 
# -------------------------
"""
    closure(u::AbstractVector, equations::GramianMomentEquations1D)

Given convective moments u[1..N] (N = M+1), compute the closure u_{N+1}
using Grad's closure.

Procedure:
    - Construct Grad's distribution function with unknowns α_k:
        f_G(c) = f_M(c; ρ,v,θ) * (1 + Σ_{k=0..M} α_k c^k)
    - Compute α by solving the linear system arising from moment matching:
        ∫ f_G(c) c^i dc = u_i, i=0..M
    - Compute closure moment:
        u_{M+1}^G = ∫ f_G(c) c^{M+1} dc
"""
function closure(u::AbstractVector, equations::GramianMomentEquations1D, ::Val{:Grad})
    M = length(u)-1                 # this is M+1 typically
    
    ρ = u[1]
    v = u[2] / ρ
    Θ = (u[3] / ρ) - v^2 # ? Check this

    # Solve for Grad coefficients α (supports Dual types)
    α = solve_alpha(u, ρ, v, Θ)

    # Compute closure: u_{M+1}^G = ∫ f_G c^{M+1} dc
    # with f_G = f_M (1 + Σ_{k=0..M} α_k c^k)
    T = promote_type(eltype(u), typeof(ρ), typeof(v), typeof(Θ))
    nquad = 2*(M+2)
    xF, wF = gausshermite(nquad)
    x = T.(xF); w = T.(wF)
    s = sqrt(T(2) * T(Θ))
    c = T(v) .+ s .* x
    pref = T(ρ) / sqrt(2 * T(pi) * T(Θ))

    return @views(
        pref * sum(w .* (c .^ (M+1)) .* (1 .+ sum(α[k+1] .* (c .^ k) for k in 0:M)))
    )
end

"""
    solve_alpha(u::AbstractVector, ρ::Real, v::Real, θ::Real; λ::Float64=0.0)

    Given convective moments u[1..N], compute the Grad coefficients α_k by solving the linear system arising from moment matching:
        ∫ f_G(c) c^i dc = u_i, i=0..M
    where f_G(c) = f_M(c; ρ,v,θ) * (1 + Σ_{k=0..M} α_k c^k)
"""
function solve_alpha(u::AbstractVector, ρ::Real, v::Real, θ::Real; λ::Float64=0.0)
    # Promote to a common numeric type (supports ForwardDiff.Dual)
    T = promote_type(eltype(u), typeof(ρ), typeof(v), typeof(θ))
    M = length(u)-1
    nquad = 2*(M+1)
    # Gauss-Hermite nodes/weights for ∫ e^{-x^2} g(x) dx
    xF, wF = gausshermite(nquad)
    x = T.(xF); w = T.(wF)
    # Construct c_j = v + sqrt(2θ)*x_j
    c = T(v) .+ sqrt(T(2.0) * T(θ)) .* x
    # Using transform with GH weights: ∫ f_M(c) c^k dc = (ρ/√π) ∑ w_j (c_j^k)
    pref = T(ρ) / sqrt(2.0 * T(pi) * T(θ))

    # Precompute powers of c up to degree 2M (needed for c^{i+k})
    maxpow = 2*M
    Cpow = ones(T, length(c), maxpow+1)
    @inbounds for p in 1:maxpow
        Cpow[:, p+1] .= Cpow[:, p] .* c
    end

    # b_i = pref * sum_j w_j * c_j^i
    b = zeros(T, M+1)
    @inbounds for i in 0:M
        b[i+1] = pref * sum(w .* Cpow[:, i+1])
    end

    # A_{ik} = pref * sum_j w_j * c_j^{i+k}
    A = zeros(T, M+1, M+1)
    @inbounds for i in 0:M
        for k in 0:M
            A[i+1, k+1] = pref * sum(w .* Cpow[:, i+k+1])
        end
    end

    # Right-hand side r = u - b (keep type T / Duals if present)
    r = u .- b

    # Regularization (Tikhonov): (A^T A + λ I) α = A^T r
    λT = T(λ)
    if λT == zero(T)
        α = A \ r
    else
        ATA = A' * A
        rhs = A' * r
        α = (ATA + λT * I) \ rhs
    end

    return α
end

# -------------------------
# Maximum entropy closure (discrete, uniform grid, exponential family)
# -------------------------

"""
    closure(u::AbstractVector, equations::GramianMomentEquations1D, ::Val{:MaxEnt})

Compute u_{M+1} by reconstructing a discrete MaxEnt distribution on a velocity grid
that matches the given convective moments u[1..M+1]. We solve for Lagrange multipliers α
in the exponential family f_i = exp(α_0 + α_1 c_i + ... + α_M c_i^M), subject to
moment constraints A f = u, then compute the next moment via the same quadrature.

Notes:
- Grid: c ∈ [v - L√Θ, v + L√Θ], uniform with n cells; defaults L=6, n=200.
- Type-generic and AD-safe (works with ForwardDiff.Dual).
"""
function closure(u::AbstractVector, equations::GramianMomentEquations1D, ::Val{:MaxEnt})
    M = length(u) - 1
    T = eltype(u)
    ρ = u[1]
    v = u[2] / ρ
    Θ = (u[3] / ρ) - v^2
    Θ = Θ > zero(T) ? Θ : T(eps(Float64))

    # Default grid parameters
    L = T(6)        # width in std devs
    n = 200         # number of cells
    cmin = v - L * sqrt(Θ)
    cmax = v + L * sqrt(Θ)

    # Solve for multipliers α matching moments 0..M
    α = _maxent_solve_multipliers(u, M, cmin, cmax, n)

    # Compute closure u_{M+1}
    dc = (cmax - cmin) / T(n)
    # Use nodes i=0..n (n+1 points) as in the Mathematica snippet
    c = [cmin + dc * T(i) for i in 0:n]
    # f_i from α
    # build polynomial α_0 + Σ α_k c^k
    poly(i) = begin
        ci = c[i]
        acc = α[1]            # α_0
        @inbounds for k in 1:M
            acc += α[k+1] * (ci^k)
        end
        acc
    end
    f = similar(c)
    @inbounds for i in eachindex(c)
        f[i] = exp(poly(i))
    end
    # next moment via rectangular rule with Dc weights
    # u_{M+1} ≈ Σ dc c_i^{M+1} f_i
    next_moment = zero(T)
    @inbounds for i in eachindex(c)
        next_moment += dc * (c[i]^(M+1)) * f[i]
    end
    return next_moment
end

# Solve for α ∈ R^{M+1} such that A(α) f(α) = u, where
# f_i(α) = exp(α_0 + α_1 c_i + ... + α_M c_i^M),
# A_{j,i} = dc * (c_i^j) with the convention c^0 ≡ 1.
function _maxent_solve_multipliers(u::AbstractVector, M::Integer, cmin, cmax, n::Integer)
    T = eltype(u)
    dc = (cmax - cmin) / T(n)
    c = [cmin + dc * T(i) for i in 0:n]   # n+1 points

    # Build A operator as a closure to avoid materializing a large matrix
    function Af(α)
        # compute f_i = exp(α⋅φ(c_i)) and moments m_j = Σ dc c_i^j f_i
        m = zeros(T, M+1)
        @inbounds for i in eachindex(c)
            # φ(c_i) polynomial
            ci = c[i]
            s = α[1]
            for k in 1:M
                s += α[k+1] * (ci^k)
            end
            fi = exp(s)
            # accumulate moments
            pow = one(T)            # ci^0
            for j in 0:M
                m[j+1] += dc * pow * fi
                pow *= ci
            end
        end
        return m
    end

    # Residual R(α) = Af(α) - u[1:M+1]
    u_vec = @views u[1:M+1]
    R(α) = Af(α) .- u_vec

    # Newton iteration
    α = zeros(T, M+1)
    maxit = 40
    rtol = T(1e-10)
    atol = T(1e-12)
    for it in 1:maxit
        r = R(α)
        rnorm = sqrt(sum(abs2, r))
        if rnorm <= rtol * max(T(1), sqrt(sum(abs2, u_vec))) + atol
            return α
        end
        # Jacobian via ForwardDiff
        J = ForwardDiff.jacobian(R, α)
        Δ = - (J \ r)
        α .= α .+ Δ
    end
    # If not converged, return the best effort α
    return α
end


# compute Gramian matrix G_n
"""
    gramian(u, n::Int)

    Gramian matrix
    G_{ij} = u_{i+j}
"""
function gramian(u, n::Int)
    G = zeros(eltype(u), n+1,n+1)
    for i=1:n+1
        for j=1:n+1
            G[i,j] = u[i+j-1]
        end
    end
    return G
end



"""
    relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}

    Two implementations, need to check which one is correct

    RHS = -1/τ (u - u_eq)
    RHS = 1/Kn (u - u_eq) (or with -?)
"""
function relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    #return SVector(ntuple(i->0.0, Mp1)...)
    # TODO: the relaxation terms need to be towards the equilibrium moments. Does the following make sense?
    prim = cons2prim(u, equations); θ = prim[3]
    eq_moments = MVector{Mp1, Float64}(undef); eq_moments[1:3] = prim[1:3]
    # compute equilibrium values for the higher moments
    val = 0.5
    for k=4:2:Mp1-1
        val = val*0.5*(k-1)
        eq_moments[k+1] = val* (2*θ)^(k/2) # gaussian integrals ∫ C^k exp(...)dC
    end
    for k=3:2:Mp1-1 eq_moments[k+1] = 0.0 end
    # τ = 1/Kn
    return SVector{Mp1}(-equations.inv_Kn .* (u .- prim2cons(eq_moments, equations)))
end

# Zero source, no RHS
zero_source(u, x, t, eqns::EqT) where {N, EqT <: Trixi.AbstractEquations{1, N}} = SVector{N}(ntuple(i->0.0, N))


# computes c^n = v^n + 0 + (n over 2) v^(n-2)C^2 + ... + (n over 1) vC^(n-1) + C^n, where u[k] is the primitive moment arising from C^(k-1) 
function binomial_moment_sum(u, n::Integer=length(u)-1)
    @assert length(u) > n > 1
    ρ = u[1]; v = u[2]
    res = v^n
    for k=2:n
        res += binomial(n, k)*v^(n-k)*u[k+1]
    end
    return res
end



# convert primitive [1, v, C^2, ...] variables to conservative [1, c, c^2, ....] variables
function moment_prim2cons(u_prim)
    m = length(u_prim)
    @assert m > 2
    cons = zeros(eltype(u_prim), m)
    # u_prim[1] is rho and u_prim[2] is v
    cons[1] = u_prim[1]; cons[2] = u_prim[1]*u_prim[2]
    for k=3:m
        cons[k] = u_prim[1]*binomial_moment_sum(u_prim, k-1)
    end
    return cons
end

# convert conservative to primitive variables
function moment_cons2prim(u_cons)
    m = length(u_cons)
    @assert m > 2
    prim = zeros(eltype(u_cons), m)
    prim[1] = u_cons[1]; prim[2] = u_cons[2]/prim[1]
    for k=3:m
        prim[k] = u_cons[k]/prim[1]- binomial_moment_sum(prim, k-1)
    end
    return prim
end

Trixi.prim2cons(u, eqns::GramianMomentEquations1D) = moment_prim2cons(u)
Trixi.cons2prim(u, eqns::GramianMomentEquations1D) = moment_cons2prim(u)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D) = u


# Derivative of closure
function dCdu(u, equations::GramianMomentEquations1D)
    # automatic differentiation
    closure_wrapped(x) = closure(x, equations)
    grad = ForwardDiff.gradient(closure_wrapped, u)
    return ForwardDiff.value(grad)
end


# Jacobian of the flux
function flux_jacobian(u, equations::GramianMomentEquations1D)
    m = length(u)
    A = zeros(eltype(u), m, m)
    for i=1:m-1 A[i, i+1] = 1 end
    A[end, :] .= dCdu(u, equations)
    
    return A
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D)
    λ_l = Trixi.max_abs_speeds(u_l, equations)
    λ_r = Trixi.max_abs_speeds(u_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end


function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D)
    # estimate the flux Jacobian eigenvalues by means of Gerschgorin
    #return max(1.0, sum(abs.(dCdu(u))))
    return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
end
