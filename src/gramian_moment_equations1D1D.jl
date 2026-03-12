
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
struct GramianMomentEquations1D{Mp1, N, RealT <: Real} <: GramianMomentEquations{Mp1, N, RealT}
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
        else
            error("Unknown closure type: $closure. Supported types are \"Gram\", \"ExtGram\", and \"Grad\".\n.")
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


"""
    relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}

    RHS = 1/Kn (u - u_eq)
"""
function relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    prim = cons2prim(u, equations); θ = prim[3]
    # Vector for equilibrium moments (primitive); the first three moments are conserved
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


"""
    computes c^n = v^n + 0 + (n over 2) v^(n-2)C^2 + ... + (n over 1) vC^(n-1) + C^n, where u[k] is the primitive moment arising from C^(k-1) 
"""
function binomial_moment_sum(u, n::Integer=length(u)-1)
    @assert length(u) > n > 1
    ρ = u[1]; v = u[2]
    res = v^n
    for k=2:n
        res += binomial(n, k)*v^(n-k)*u[k+1]
    end
    return res
end



"""
    convert primitive [1, v, C^2, ...] variables to conservative [1, c, c^2, ....] variables
"""
function moment_prim2cons(u_prim, eqns::GramianMomentEquations1D)
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


"""
    convert conservative to primitive variables
"""
function moment_cons2prim(u_cons, eqns::GramianMomentEquations1D)
    m = length(u_cons)
    @assert m > 2
    prim = zeros(eltype(u_cons), m)
    prim[1] = u_cons[1]; prim[2] = u_cons[2]/prim[1]
    for k=3:m
        prim[k] = u_cons[k]/prim[1]- binomial_moment_sum(prim, k-1)
    end
    return prim
end

Trixi.prim2cons(u, eqns::GramianMomentEquations1D) = moment_prim2cons(u, eqns)
Trixi.cons2prim(u, eqns::GramianMomentEquations1D) = moment_cons2prim(u, eqns)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D) = u


"""
    Derivative of the closure function w.r.t. the moments u
"""
function dCdu(u, equations::GramianMomentEquations1D)
    # automatic differentiation
    closure_wrapped(x) = closure(x, equations)
    grad = ForwardDiff.gradient(closure_wrapped, u)
    return ForwardDiff.value(grad)
end


"""
    Jacobian of the flux function
"""
function flux_jacobian(u, equations::GramianMomentEquations1D)
    m = length(u)
    A = zeros(eltype(u), m, m)
    for i=1:m-1 A[i, i+1] = 1 end
    A[end, :] .= dCdu(u, equations)
    
    return A
end

"""
    Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
"""
function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D)
    λ_l = Trixi.max_abs_speeds(u_l, equations)
    λ_r = Trixi.max_abs_speeds(u_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end


"""
    Estimate the flux Jacobian eigenvalues by means of Gerschgorin
"""
function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D)
    return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
end
