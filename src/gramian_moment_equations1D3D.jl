
# (note: replacing all zeroes() calls with corresponding SVector and MVector calls would speed up things at the cost of arbitrary vector input types)

# * Basic struct for equations
"""
struct GramianMomentEquations1D3D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}

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
        - for the moment hard-coded to M=4
- `Knudsen::Real`: Knudsen number for the system
- `closure::Union{Symbol,String}`: Closure type, can be "Gram", "ExtGram", or "Grad"
# Notes
- The system has `M+1` variables in total
- The parameter `χ` is calculated differently depending on whether `M` is even or odd
"""
struct GramianMomentEquations1D3D{Mp1, N, RealT <: Real} <: Trixi.AbstractEquations{1, Mp1}
    inv_Kn::RealT # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    n::Int
    closure::Symbol

    function GramianMomentEquations1D3D(M::Integer, Knudsen::Real, closure::String="ExtGram"; χ_set="optimal")
        @assert M > 1
        @assert M == 4 "Currently only M=4 is supported."
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
        else
            error("Unknown closure type: $closure. Supported types are \"Gram\" and \"ExtGram\".\n.")
        end
        new{35, n, typeof(Knudsen)}(inv(Knudsen), χ, n, closure_value)
    end
end

# * Conservative to Primitive variables and vice versa
# ToDo: Change for 3D velocity case
Trixi.varnames(::typeof(cons2cons), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = ntuple(i->"u^($(i-1))", Mp1)
Trixi.varnames(::typeof(cons2prim), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = ntuple(i->"w^($(i-1))", Mp1)

# * basic variable definitions to be used as e.g. sc indicator variable
Trixi.density(u, eqns::GramianMomentEquations1D3D{Mp1}) where {Mp1} = u[1]


"""
    Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
    Flux function: F(U) = u_{k+1}, k=0,…,M where u_{M+1} = C(u0, u1, …, uM)
"""
function Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
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
closure(u, equations::GramianMomentEquations1D3D) = closure(u, equations, Val(equations.closure))

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
    relaxation_source(u, x, t, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}

    RHS = 1/Kn (u - u_eq) (or with -?)
"""
# todo: implement equilibrium moments for 3D velocity case
# function relaxation_source(u, x, t, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
#     #return SVector(ntuple(i->0.0, Mp1)...)
#     # TODO: the relaxation terms need to be towards the equilibrium moments. Does the following make sense?
#     prim = cons2prim(u, equations); θ = prim[3]
#     eq_moments = MVector{Mp1, Float64}(undef); eq_moments[1:3] = prim[1:3]
#     # compute equilibrium values for the higher moments
#     val = 0.5
#     for k=4:2:Mp1-1
#         val = val*0.5*(k-1)
#         eq_moments[k+1] = val* (2*θ)^(k/2) # gaussian integrals ∫ C^k exp(...)dC
#     end
#     for k=3:2:Mp1-1 eq_moments[k+1] = 0.0 end
#     # τ = 1/Kn
#     return SVector{Mp1}(-equations.inv_Kn .* (u .- prim2cons(eq_moments, equations)))
# end

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

# ToDo: Need to check for 3D velocity case
Trixi.prim2cons(u, eqns::GramianMomentEquations1D3D) = moment_prim2cons(u)
Trixi.cons2prim(u, eqns::GramianMomentEquations1D3D) = moment_cons2prim(u)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D3D) = u


# Derivative of closure
function dCdu(u, equations::GramianMomentEquations1D3D)
    # automatic differentiation
    closure_wrapped(x) = closure(x, equations)
    grad = ForwardDiff.gradient(closure_wrapped, u)
    return ForwardDiff.value(grad)
end


# Jacobian of the flux
function flux_jacobian(u, equations::GramianMomentEquations1D3D)
    m = length(u)
    A = zeros(eltype(u), m, m)
    for i=1:m-1 A[i, i+1] = 1 end
    A[end, :] .= dCdu(u, equations)
    
    return A
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3D)
    λ_l = Trixi.max_abs_speeds(u_l, equations)
    λ_r = Trixi.max_abs_speeds(u_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end


function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D3D)
    # estimate the flux Jacobian eigenvalues by means of Gerschgorin
    #return max(1.0, sum(abs.(dCdu(u))))
    return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
end
