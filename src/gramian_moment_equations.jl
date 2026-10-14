
using StaticArrays: SUnitRange
"""
This file defines the GramianMomentEquations type and the 1D-1D closure relations
The GramianMomentEquations type is a subtype of Trixi.AbstractEquations and represents the system of equations for the Gramian moments.
"""


"""
Abstract type for the Gramian Moment Equations. This type is a subtype of Trixi.AbstractEquations and represents the system of equations for the Gramian moments.

This class serves as a marker for the Gramian Moment Equations and does not contain any fields. The specific implementation of the equations is done in the subtypes, such as GramianMomentEquations and GramianMomentEquations3D.
It allows to use different closures.
"""
abstract type GramianMomentEquations{Mp1, N, RealT <: Real, C} <: Trixi.AbstractEquations{1, Mp1} end

"""
    closure(u, equations::GramianMomentEquations)

    Decides which closure implementation is used based on the closure type parameter `C` of the equations
"""
@inline function closure(u, equations::GramianMomentEquations{Mp1, N, RealT, C}) where {Mp1, N, RealT, C}
    return closure(u, equations, Val(C))
end

# -------------------------
# Gramian closure 
# -------------------------

# Even case (classical and extended)
"""
    closure(u, equations::GramianMomentEquations, ::Val{:GramEven})

    Gramian closure for the even case
"""
function closure(u, equations::GramianMomentEquations{Mp1, n}, ::Val{:GramEven}) where {Mp1, n}
    @assert iseven(length(u) - 1) # u[0, ..., M]
    invG_nm1 = inv(gramian(u, Val(n-1)))
    return u[SUnitRange(n+2, 2n+1)]' * invG_nm1 * u[SUnitRange(n+1, 2n)]
end

"""
    closure(u, equations::GramianMomentEquations, ::Val{:ExtGramEven})

    Extended Gramian closure for the even case
"""
function closure(u, equations::GramianMomentEquations{Mp1, n}, ::Val{:ExtGramEven}) where {Mp1, n}
    @assert iseven(length(u) - 1) # u[0, ..., M]
    invG_nm1 = inv(gramian(u, Val(n-1)))
    invG_nm2 = inv(gramian(u, Val(n-2)))
    u_np_2n    = u[SUnitRange(n+1, 2n)]   # (u_n, …, u_{2n-1})
    u_np2_2np1  = u[SUnitRange(n+2, 2n+1)] # (u_{n+1}, …, u_{2n})
    u_n_2nm2  = u[SUnitRange(n, 2n-2)]   # (u_{n-1}, …, u_{2n-3})
    u_np1_2nm1 = u[SUnitRange(n+1, 2n-1)] # (u_n, …, u_{2n-2})
    return (
        u_np2_2np1' * invG_nm1 * u_np_2n
    ) + equations.χ * (
        u[2n+1] - u_np_2n' * invG_nm1 * u_np_2n
    ) / (
        u[2n-1] - u_n_2nm2' * invG_nm2 * u_n_2nm2
    ) * (
        u[2n] - u_np1_2nm1' * invG_nm2 * u_n_2nm2
    )
end

# Odd case (classical and extended)
"""
    closure(u, equations::GramianMomentEquations, ::Val{:GramOdd})

    Gramian closure for the odd case
"""
function closure(u, equations::GramianMomentEquations{Mp1, n}, ::Val{:GramOdd}) where {Mp1, n}
    @assert isodd(length(u) - 1) # u[0, ..., M]
    invG_nm1 = inv(gramian(u, Val(n-1)))
    u_n = u[SUnitRange(n+1, 2n)] # (u_n, …, u_{2n-1})
    return u_n' * invG_nm1 * u_n
end

"""
    closure(u, equations::GramianMomentEquations, ::Val{:GramOdd})

    Extended Gramian closure for the odd case
"""
sigma_nn(u, ::Val{n}) where {n}   = u[2n+1] - u[SUnitRange(n+1, 2n)]' * inv(gramian(u, Val(n-1))) * u[SUnitRange(n+1, 2n)]
sigma_nmnm(u, ::Val{n}) where {n} = u[2n-1] - u[SUnitRange(n, 2n-2)]' * inv(gramian(u, Val(n-2))) * u[SUnitRange(n, 2n-2)]
sigma_nmnm(u, ::Val{1}) = u[1] # empty Gramian G_{-1}: the Schur complement is just u_0
sigma_nnp(u, ::Val{n}) where {n}  = u[2n+2] - u[SUnitRange(n+2, 2n+1)]' * inv(gramian(u, Val(n-1))) * u[SUnitRange(n+1, 2n)]
sigma_nmn(u, ::Val{n}) where {n}  = u[2n] - u[SUnitRange(n+1, 2n-1)]' * inv(gramian(u, Val(n-2))) * u[SUnitRange(n, 2n-2)]

# Recurrence coefficients α_n, β_n of the orthogonal polynomials
function an(u, ::Val{n}) where {n}
    if n == 0
        return u[2] / u[1]
    elseif n == 1
        return sigma_nnp(u, Val(1)) / sigma_nn(u, Val(1)) - u[2] / u[1]
    else
        return sigma_nnp(u, Val(n)) / sigma_nn(u, Val(n)) - sigma_nmn(u, Val(n)) / sigma_nmnm(u, Val(n))
    end
end
bn(u, ::Val{0}) = zero(eltype(u))
bn(u, ::Val{n}) where {n} = sigma_nn(u, Val(n)) / sigma_nmnm(u, Val(n))

# (α_0, …, α_{n-1}) and (β_0, …, β_{n-1}) as tuples. Built recursively so that every `Val(k)` is a
# compile-time constant (an `ntuple` over `k -> an(u, Val(k-1))` would dispatch at run time and allocate).
recurrence_coefficients(u, ::Val{0}) = ((), ())
function recurrence_coefficients(u, ::Val{n}) where {n}
    α, β = recurrence_coefficients(u, Val(n-1))
    return ((α..., an(u, Val(n-1))), (β..., bn(u, Val(n-1))))
end

function closure(u, equations::GramianMomentEquations{Mp1, n}, ::Val{:ExtGramOdd}) where {Mp1, n}
    @assert isodd(length(u) - 1) # u[0, ..., M]

    # α_k, β_k for k = 0, …, n-1
    α, β = recurrence_coefficients(u, Val(n))

    β_k = zero(eltype(u))
    for k in 1:n-1
        β_k += 2 / (n-1) * β[k+1]
    end
    for l in 0:n-1
        β_k += 1 / (n-1) * (α[n] - α[l+1])^2
    end

    u_n = u[SUnitRange(n+1, 2n)] # (u_n, …, u_{2n-1})
    return u_n' * inv(gramian(u, Val(n-1))) * u_n + sigma_nmnm(u, Val(n)) * β_k
end


# -------------------------
# Grad closure 
# -------------------------
"""
    closure(u::AbstractVector, equations::GramianMomentEquations)

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
function closure(u::AbstractVector, equations::GramianMomentEquations, ::Val{:Grad})
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



# compute Gramian matrix G_n
"""
    gramian(u, n::Int)

    Gramian matrix
    G_{ij} = u_{i+j}
"""
@inline function gramian(u, ::Val{n}) where {n}
    G = MMatrix{n+1, n+1, eltype(u)}(undef)
    for j in 1:n+1, i in 1:n+1
        G[i, j] = u[i+j-1]
    end
    return SMatrix(G)
end