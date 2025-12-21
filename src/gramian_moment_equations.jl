
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
    """
    First MP1-1 flux components from shifed moments
    ∂_t u^k + ∂_x u^{k+1} = ... (0)
    Last u from closure
    """
    return SVector(ntuple(i->u[i+1], Mp1-1)..., closure(u, equations, Val(equations.closure)))
end


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
# -------------------------
# closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGram}) = closure(u, equations, Val(iseven(length(u)-1) ? :ExtGramEven : :ExtGramOdd))

# Even case (classical and extended)
function closure(u, equations::GramianMomentEquations1D, ::Val{:GramEven})
    """
    closure(u, equations::GramianMomentEquations1D, ::Val{:GramEven})

    Gramian closure for the even case
    """
    M = length(u)-1 # u[0, ..., M]
    @assert iseven(M)
    n = equations.n 
    invG_nm1 = inv(gramian(u, n-1))
    return @views(
        u[n+2:2n+1]'*invG_nm1*u[n+1:2n]
    )
end

function closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGramEven})
    """
    closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGramEven})

    Extended Gramian closure for the even case
    """
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
function closure(u, equations::GramianMomentEquations1D, ::Val{:GramOdd})
    """
    closure(u, equations::GramianMomentEquations1D, ::Val{:GramOdd})

    Gramian closure for the odd case
    """
    M = length(u)-1 # u[0, ..., M]
    @assert isodd(M)
    n = equations.n
    invG_nm1 = inv(gramian(u, n-1))
    return @views (
        u[n+1:2n]'*invG_nm1*u[n+1:2n]
    )
end

function closure(u, equations::GramianMomentEquations1D, ::Val{:ExtGramOdd})
    """
    closure(u, equations::GramianMomentEquations1D, ::Val{:GramOdd})

    Extended Gramian closure for the odd case
    """
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
function closure(u::AbstractVector, equations::GramianMomentEquations1D, ::Val{:Grad})
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

function solve_alpha(u::AbstractVector, ρ::Real, v::Real, θ::Real; λ::Float64=0.0)
    """
        solve_alpha(u::AbstractVector, ρ::Real, v::Real, θ::Real; λ::Float64=0.0)

        Given convective moments u[1..N], compute the Grad coefficients α_k by solving the linear system arising from moment matching:
            ∫ f_G(c) c^i dc = u_i, i=0..M
        where f_G(c) = f_M(c; ρ,v,θ) * (1 + Σ_{k=0..M} α_k c^k)
    """
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
function gramian(u, n::Int)
    """
        gramian(u, n::Int)

        Gramian matrix
        G_{ij} = u_{i+j}
    """
    G = zeros(eltype(u), n+1,n+1)
    for i=1:n+1
        for j=1:n+1
            G[i,j] = u[i+j-1]
        end
    end
    return G
end



function relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    """
        relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}

        RHS = 1/Kn (u - u_eq)
    """
    prim = cons2prim(u, equations); θ = prim[3]
    eq_moments = MVector{Mp1, Float64}(undef); eq_moments[1:3] = prim[1:3]
    # compute equilibrium values for the higher moments
    val = 0.5
    for k=4:2:Mp1-1
        val = val*0.5*(k-1)
        eq_moments[k+1] = val* (2*θ)^(k/2) # gaussian integrals ∫ C^k exp(...)dC
    end
    for k=3:2:Mp1-1 eq_moments[k+1] = 0.0 end
    return SVector{Mp1}(-equations.inv_Kn .* (u .- prim2cons(eq_moments, equations)))
end

""" 
    Zero source, no RHS
"""
zero_source(u, x, t, eqns::EqT) where {N, EqT <: Trixi.AbstractEquations{1, N}} = SVector{N}(ntuple(i->0.0, N))


function binomial_moment_sum(u, n::Integer=length(u)-1)
    """
        computes c^n = v^n + 0 + (n over 2) v^(n-2)C^2 + ... + (n over 1) vC^(n-1) + C^n, where u[k] is the primitive moment arising from C^(k-1) 
    """
    @assert length(u) > n > 1
    ρ = u[1]; v = u[2]
    res = v^n
    for k=2:n
        res += binomial(n, k)*v^(n-k)*u[k+1]
    end
    return res
end



function moment_prim2cons(u_prim)
    """
        convert primitive [1, v, C^2, ...] variables to conservative [1, c, c^2, ....] variables
    """
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


function moment_cons2prim(u_cons)
    """
        convert conservative to primitive variables
    """
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


function dCdu(u, equations::GramianMomentEquations1D)
    """
        Derivative of the closure function w.r.t. the moments u
    """
    # automatic differentiation
    closure_wrapped(x) = closure(x, equations)
    grad = ForwardDiff.gradient(closure_wrapped, u)
    return ForwardDiff.value(grad)
end


function flux_jacobian(u, equations::GramianMomentEquations1D)
    """
        Jacobian of the flux function
    """
    m = length(u)
    A = zeros(eltype(u), m, m)
    for i=1:m-1 A[i, i+1] = 1 end
    A[end, :] .= dCdu(u, equations)
    
    return A
end

function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D)
    """
        Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
    """
    λ_l = Trixi.max_abs_speeds(u_l, equations)
    λ_r = Trixi.max_abs_speeds(u_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end


function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D)
    """
        Estimate the flux Jacobian eigenvalues by means of Gerschgorin
    """
    return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
end
