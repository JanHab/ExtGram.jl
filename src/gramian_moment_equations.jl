
# (note: replacing all zeroes() calls with corresponding SVector and MVector calls would speed up things at the cost of arbitrary vector input types)

# * Basic struct for equations
"""
GramianMomentEquations1D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}

A struct representing one-dimensional Gramian moment equations for kinetic theory.

# Fields
- `inv_Kn::RealT`: The inverse of the Knudsen number, used in the relaxation source term
- `χ::RealT`: A scaling parameter calculated based on the number of moments
- `n::Int`: Internal parameter derived from the number of moments
- `extended::Bool`: Flag indicating whether to use extended formulation

# Constructor
    GramianMomentEquations1D(M::Integer, Knudsen::Real, extended=true::Bool)

Constructs a `GramianMomentEquations1D` system with `M` moments and specified Knudsen number.

# Arguments
- `M::Integer`: Number of moments (must be greater than 1)
- `Knudsen::Real`: Knudsen number for the system
- `extended::Bool=true`: Whether to use the extended formulation

# Notes
- The system has `M+1` variables in total
- The parameter `χ` is calculated differently depending on whether `M` is even or odd
"""
struct GramianMomentEquations1D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}
    inv_Kn::RealT # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    n::Int
    extended::Bool

    function GramianMomentEquations1D(M::Integer, Knudsen::Real, extended=true::Bool; χ_set="optimal")
        @assert M > 1 # todo: remove this later
        if iseven(M)
            n = Int(M/2)
            χ = (n+1)/n
        else
            n = Int((M+1)/2)
            χ = (n+1)/(2n)
        end
        # overwrite if χ_set is given explicitly
        if χ_set !== "optimal"
            χ = χ_set
        end
        new{M+1, n, typeof(Knudsen)}(inv(Knudsen), χ, n, extended)
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
    return SVector(ntuple(i->u[i+1], Mp1-1)..., closure(u, equations; verbose_=true))
end
isinvertible(A::Matrix{Float64}) = !isapprox(det(BigFloat.(A)), 0, atol = 1e-18)


########################## Closure ##########################

# Overloading to determine even/odd case
"""
    closure(u, equations::GramianMomentEquations1D; verbose_::Bool=false)

    Decides which closure implementation is used based on the number of Moments (even vs. odd)
"""
closure(u, equations::GramianMomentEquations1D; verbose_::Bool=false) = grad_closure_convective(u, equations)


# # -------------------------
# # Build Atmp via Gauss-Hermite quadrature (no sympy)
# # -------------------------
# """
#     build_Atmp(Mmax; gh_n=2*(Mmax+2))

# Build the (Mmax+2) x (Mmax+2) matrix Atmp with entries
#     Atmp[i,j] = ∫_{-∞}^{∞} (1/√(2π)) e^{-ξ^2/2} ξ^{(i+j-2)} dξ
# computed by Gauss-Hermite quadrature.

# We use the substitution ξ = sqrt(2) ζ so that
#     ∫ e^{-ξ^2/2} g(ξ) dξ = sqrt(2) ∫ e^{-ζ^2} g(√2 ζ) dζ
# and FastGaussQuadrature.gausshermite returns nodes ζ_k and weights w_k for ∫ e^{-ζ^2} h(ζ) dζ ≈ Σ w_k h(ζ_k).
# """
# function build_Atmp(Mmax)
#     # nodes ξ and weights w approximate ∫_{-∞}^{∞} e^{-ξ^2} h(ξ) dξ
#     # ξ, w = gausshermite(Mmax+1) # +1 for good measure, should not be necessary
#     # sizeA = Mmax + 2
#     # Atmp = zeros(Float64, sizeA, sizeA)

#     # # # p = i+j-2 is the exponent for ξ^p
#     # # # A_{ij} = (1/√(2π)) * 2^{(p+1)/2} * ∫ e^{-ξ^2} ξ^p dξ
#     # # for i in 1:sizeA
#     # #     for j in 1:sizeA
#     # #         p = i + j - 2
#     # #         # approximate integral
#     # #         S = sum(w .* (ξ .^ p))
#     # #         Atmp[i, j] = 1 / sqrt(2π) * 2^((p+1)/2) * S
#     # #     end
#     # # end
#     # # vectorized version
#     # fg = 1 / sqrt(2π) * exp.(- ξ.^2 / 2) .* [sum(w .* (ξ .^ k)) for k in 0:(Mmax+1)] # !0:(2*sizeA-2)] # factor outside integral
#     # for i in 1:sizeA
#     #     for j in 1:sizeA
#     #         p = i + j - 2
#     #         Atmp[i, j] = fg[p+1] * 2^((p+1)/2)
#     #     end
#     # end
#     # return Atmp

#     # nodes and weights for ∫ e^{-x^2} f(x) dx  (physicists' normalization)
#     ξ, w = gausshermite(2*Mmax + 4)  # extra accuracy

#     # we need ∫ e^{-ξ²/2} ξ^p ξ^q dξ
#     # but gausshermite() integrates e^{-ξ²} f(ξ), so change variable:
#     # ∫ e^{-ξ²/2} f(ξ) dξ = ∫ e^{-ξ²} f(ξ) * e^{ξ²/2} dξ
#     # weightfix = exp.(ξ.^2 ./ 2) ./ sqrt(2π)
#     weightfix = exp.(ξ.^2 ./ 2) ./ sqrt(2π) # ? Which sign do we need here?

#     # Build A matrix: A[i,j] = ∫ ξ^(i-1) * ξ^(j-1) * e^{-ξ²/2}/√(2π) dξ
#     # = expectation of ξ^(i+j-2)
#     A = zeros(Float64, Mmax + 2, Mmax + 2)
#     for i in 1:Mmax+2, j in 1:Mmax+2
#         A[i,j] = sum(w .* weightfix .* (ξ .^ (i + j - 2)))
#     end

#     return A
# end

# # -------------------------
# # Compute NextGrad via linear algebra
# # -------------------------
# """
#     compute_NextGrad(Mmax; gh_n=...)

# Return a Dict mapping M -> NextGrad vector for M = 4..Mmax.
# Each NextGrad[M] is a 1×(M+1) row vector (as Vector{Float64}) such that
#   u_{M+1} = NextGrad[M] * u_{0:M}
# in the dimensionless central coordinate system used by Grad.
# """
# function compute_NextGrad(Mmax)
#     Atmp = build_Atmp(Mmax)
#     NextGrad = Dict{Int, Vector{Float64}}()
#     for M in 4:Mmax
#         A11 = Atmp[1:(M+1), 1:(M+1)]
#         A21 = Atmp[M+2, 1:(M+1)]          # row vector
#         # small sizes: computing inv is fine; optionally do \ for stability
#         coeffs = A21' * inv(A11) # ? Do we take the inverse? Is mathematica using rows first and julia columns first or other way around?
#         NextGrad[M] = vec(coeffs)
#     end
#     return NextGrad
# end

# function compute_NextGrad(N::Int)
#     # Moment function for standard normal
#     μ(n) = isodd(n) ? 0.0 : factorial(big(2*(n ÷ 2))) / (2.0^(n ÷ 2) * factorial(big(n ÷ 2)))

#     # Precompute all needed moments up to 2N+1
#     μvals = [Float64(μ(k)) for k in 0:(2N+1)]

#     # Build moment matrix A and vector b
#     A = [μvals[i+j+1] for i in 0:N, j in 0:N]  # Julia is 1-indexed
#     b = [μvals[N+1+j+1] for j in 0:N]

#     # Solve for NextGrad
#     NextGrad = b' * inv(A)
#     return vec(NextGrad)
# end

# # -------------------------
# # Grad closure function
# # -------------------------
# """
#     grad_closure_convective(m::AbstractVector, NextGrad::Dict)

# Given convective moments m[1..N] (N = M+1), compute the closure m_{N+1}
# using Grad's closure built into NextGrad.

# Procedure:
#     - compute central moments via m2rho
#     - set central first moment to zero (centering)
#     - nondimensionalize central moments: u_k = rho_k / (rho0 * Theta^{k/2})
#     - apply NextGrad[M] to u_0..u_M to get u_{M+1}
#     - dimensionalize rho_{M+1} = rho0 * Theta^{(M+1)/2} * u_{M+1}
#     - restore a and convert central->convective using rho2m and return last entry
# """
# function grad_closure_convective(u::AbstractVector, equations::GramianMomentEquations1D) #!NextGrad::Dict) # compute NextGrad up to M=length(u)-1
#     N = length(u)                 # this is M+1 typically
#     M = N - 1
#     NextGrad = compute_NextGrad(N)
#     # @assert haskey(NextGrad, M) "NextGrad for M=$M not present"

#     # convective -> central
#     w = moment_cons2prim(u)
#     ρ = w[1]; θ=w[3]/ρ

#     # # nondimensionalize: u_k = rho_k / (rho0 * Θ^{k/2}), k=0..M
#     scale = ρ .* [θ^(k/2) for k in 0:(M)]
#     w_scaled = w ./ scale

#     # compute next dimensionless central moment
#     # println("Using Grad closure for M=$M")
#     # println("NextGrad = ", NextGrad, "\n")
#     # println("Size of input moments: ", length(u), "\n")
#     # println("Input moments: ", u, "\n")
#     # println("Scaled moments = ", w_scaled, "\n")
#     # ? w_next = dot(NextGrad[M], w_scaled)   # scalar
#     w_next = dot(NextGrad[1:M+1], w_scaled)   # scalar

#     # dimensionalize: rho_{M+1} = rho0 * Θ^{(M+1)/2} * u_next
#     # append to convective moments
#     w_extended = vcat(w, ρ * θ^(N/2) * w_next)

#     # convective moments back
#     u_extended = moment_prim2cons(w_extended)

#     return u_extended[end]   # the closure u_{M+1}
# end

# """
#     maxwellian_convective_moments(ρ, v, θ, K; nquad=…)

# Compute m_k = ∫ f_M(c; ρ,v,θ) c^k dc for k=0..K using Gauss–Hermite.
# Returns a Vector{T} with length K+1, where T is a common type of inputs.
# """
# function maxwellian_convective_moments(ρ, v, θ, K; nquad::Int=0)
#     T = promote_type(typeof(ρ), typeof(v), typeof(θ))
#     nquad = nquad == 0 ? max(2*(K+1), 50) : nquad
#     xF, wF = gausshermite(nquad)             # Float64 nodes/weights
#     x = T.(xF); w = T.(wF)
#     s = sqrt(T(2) * θ)
#     c = v .+ s .* x                          # Vector{T}
#     pref = ρ / sqrt(T(pi))                   # prefactor ρ/√π

#     # build powers up to K
#     Cpow = ones(T, length(c), K+1)
#     @inbounds for p in 1:K
#         Cpow[:, p+1] .= Cpow[:, p] .* c
#     end

#     m = zeros(T, K+1)
#     @inbounds for i in 0:K
#         m[i+1] = pref * sum(w .* Cpow[:, i+1])
#     end
#     return m
# end

# function solve_alpha(u::AbstractVector, rho, v, theta; nquad::Int = 0, λ = 0)
#     # Element type to support ForwardDiff.Dual
#     T = promote_type(eltype(u), typeof(rho), typeof(v), typeof(theta))
#     M = length(u) - 1
#     nquad = nquad == 0 ? max(2*(M+1), 50) : nquad

#     # Maxwellian moments up to 2M
#     m = maxwellian_convective_moments(rho, v, theta, 2M; nquad=nquad)

#     # Build b and A with element type T
#     b = zeros(T, M+1)
#     @inbounds for i in 0:M
#         b[i+1] = m[i+1]
#     end

#     A = zeros(T, M+1, M+1)
#     @inbounds for i in 0:M
#         for k in 0:M
#             A[i+1, k+1] = m[i+k+1]
#         end
#     end

#     # Right-hand side r = u - b (no Float64 conversion; keep Duals if present)
#     r = u .- b

#     # Regularization (Tikhonov): (A^T A + λ I) α = A^T r
#     if λ == zero(T)
#         α = A \ r
#     else
#         ATA = A' * A
#         rhs = A' * r
#         α = (ATA + (λ * I)) \ rhs
#     end

#     return α
# end

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

function grad_closure_convective(u::AbstractVector, equations::GramianMomentEquations1D) #!NextGrad::Dict) # compute NextGrad up to M=length(u)-1
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

    # # Build polynomial factor (1 + Σ α_k c^k)
    # poly = ones(T, length(c))
    # @inbounds for k in 0:M
    #     poly .+= α[k+1] .* (c .^ k)
    # end

    # return pref * sum(w .* (c .^ (M+1)) .* poly)
end


# Even case
"""
    closure(u, equations::GramianMomentEquations1D, ::Val{true}; verbose_=false)

    Closure for the even case
"""
function closure(u, equations::GramianMomentEquations1D, ::Val{true}; verbose_=false)
    M = length(u)-1 # u[0, ..., M]
    @assert iseven(M)
    n = equations.n 
    invG_nm1 = inv(gramian(u, n-1))
    if equations.extended==true
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
    else
        return @views(
            u[n+2:2n+1]'*invG_nm1*u[n+1:2n]
        )
    end
end

# Odd case
"""
    closure(u, equations::GramianMomentEquations1D, ::Val{false}; verbose_=false)

    Closure for the odd case
"""
function closure(u, equations::GramianMomentEquations1D, ::Val{false}; verbose_=false)
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    if equations.extended==true
        # invG_nm2 = inv(gramian(u, n-2))
        # # A = u[n+2:2n]'*invG_nm2*u[n:2n-2]
        # # B = u[2n] - u[n+1:2n-1]'*invG_nm2*u[n:2n-2]
        # # D = u[2n-1] - u[n:2n-2]'*invG_nm2*u[n:2n-2]
        # # closure = A + equations.χ * B^2 / D
        # invG_nm1 = inv(gramian(u, n-1))
        
        # some ad-hoc odd-order closures from Morin paper
        if n==2 # M=3
            invG_1 = inv(gramian(u, 1))
            invG_0 = 1 / u[1]#!inv(gramian(u, 0))
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
            # New closure from SymPy (u6 computation with proper index shift)
            # SymPy: A = 1 + (u4 - [u2,u3]'*invG_1*[u2,u3])/(u2-u1*invG_0*u1) + 0.5*((...)^2 + (...)^2)
            # Julia: A = 1 + (u5 - [u3,u4]'*invG_1*[u3,u4])/(u3-u2*invG_0*u2) + 0.5*((...)^2 + (...)^2)
            invG_1 = inv(gramian(u, 1))
            invG_0 = 1 / u[1]#!inv(gramian(u, 0))
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
            
            # SymPy: u6 = [u3,u4,u5]'*invG_2*[u3,u4,u5] + (u4-[u2,u3]'*invG_1*[u2,u3])*A
            # Julia: u6 = [u4,u5,u6]'*invG_2*[u4,u5,u6] + (u5-[u3,u4]'*invG_1*[u3,u4])*A
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
    else
        invG_nm1 = inv(gramian(u, n-1))
        return @views (
            u[n+1:2n]'*invG_nm1*u[n+1:2n]
        )
    end
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
