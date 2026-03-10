
"""
This file defines the GramianMomentEquations type and the 1D-1D closure relations
The GramianMomentEquations type is a subtype of Trixi.AbstractEquations and represents the system of equations for the Gramian moments.
"""


"""
Abstract type for the Gramian Moment Equations. This type is a subtype of Trixi.AbstractEquations and represents the system of equations for the Gramian moments.

This class serves as a marker for the Gramian Moment Equations and does not contain any fields. The specific implementation of the equations is done in the subtypes, such as GramianMomentEquations and GramianMomentEquations3D.
It allows to use different closures.
"""
abstract type GramianMomentEquations{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}end




########################## Closure ##########################

# Overloading to determine Gram/ExtGram/Grad and even/odd case
"""
    closure(u, equations::GramianMomentEquations)

    Decides which closure implementation is used based on the equations.closure value
"""
closure(u, equations::GramianMomentEquations) = closure(u, equations, Val(equations.closure))

# -------------------------
# Gramian closure 
# -------------------------

# Even case (classical and extended)
"""
    closure(u, equations::GramianMomentEquations, ::Val{:GramEven})

    Gramian closure for the even case
"""
function closure(u, equations::GramianMomentEquations, ::Val{:GramEven})
    M = length(u)-1 # u[0, ..., M]
    @assert iseven(M)
    n = equations.n 
    invG_nm1 = inv(gramian(u, n-1))
    return @views(
        u[n+2:2n+1]'*invG_nm1*u[n+1:2n]
    )
end

"""
    closure(u, equations::GramianMomentEquations, ::Val{:ExtGramEven})

    Extended Gramian closure for the even case
"""
function closure(u, equations::GramianMomentEquations, ::Val{:ExtGramEven})
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
    closure(u, equations::GramianMomentEquations, ::Val{:GramOdd})

    Gramian closure for the odd case
"""
function closure(u, equations::GramianMomentEquations, ::Val{:GramOdd})
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    invG_nm1 = inv(gramian(u, n-1))
    return @views (
        u[n+1:2n]'*invG_nm1*u[n+1:2n]
    )
end

"""
    closure(u, equations::GramianMomentEquations, ::Val{:GramOdd})

    Extended Gramian closure for the odd case
"""
function closure(u, equations::GramianMomentEquations, ::Val{:ExtGramOdd})
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    s = 1 # shift for index as julia is 1-based

    sigma_nn(u, n) = u[2n+s] - u[n+s:2n-1+s]'*inv(gramian(u, n-1))*u[n+s:2n-1+s]
    sigma_nmnm(u, n) = u[2n-2+s] - u[n-1+s:2n-3+s]'*inv(gramian(u, n-2))*u[n-1+s:2n-3+s]
    sigma_nnp(u, n) = u[2n+1+s] - u[n+1+s:2n+s]'*inv(gramian(u, n-1))*u[n+s:2n-1+s]
    sigma_nmn(u, n) = u[2n-1+s] - u[n+s:2n-2+s]'*inv(gramian(u, n-2))*u[n-1+s:2n-3+s]

    function an(u, n)
        if n == 0
            return u[1+s] / u[0+s]
        elseif n == 1
            return sigma_nnp(u, n) / sigma_nn(u, n) - u[1+s] / u[0+s]
        else
            return sigma_nnp(u, n) / sigma_nn(u, n) - sigma_nmn(u, n) / sigma_nmnm(u, n)
        end
    end
    function bn(u, n)
        if n == 0
            return 0
        else
            return sigma_nn(u, n) / sigma_nmnm(u, n)
        end
    end

    α = zeros(eltype(u), n)
    β = zeros(eltype(u), n)
    α[0+s] = u[1+s] / u[0+s]
    β[0+s] = 0.0
    for k=1:n-1
        α[k+s] = an(u, k)
        β[k+s] = bn(u, k)
    end
    
    β_k = 0.0
    for k=1:n-1
        β_k += 2 / (n-1) * β[k+s]
    end
    for l=0:n-1
        β_k += 1 / (n-1) * (α[n-1+s] - α[l+s])^2
    end
    closure = @views(
        u[n+s:2n-1+s]' * inv(gramian(u, n-1)) * u[n+s:2n-1+s]
    ) + (
        sigma_nmnm(u, n)
    ) * (
        β_k
    )

    return closure
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
function gramian(u, n::Int)
    G = zeros(eltype(u), n+1,n+1)
    for i=1:n+1
        for j=1:n+1
            G[i,j] = u[i+j-1]
        end
    end
    return G
end


