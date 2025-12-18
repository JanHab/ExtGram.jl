"""
    Implementation of the one-dimensional Gramian moment equations with 3D velocity dependence.
"""

struct Constants{F<:Function}
    """
        Struct to hold precomputed constants for Gramian moment equations in 1D with 3D velocity.

        # Fields
        - `fp_func::F`: Compiled function for evaluating basis functions at given angles.
        - `Angles::Vector{Tuple{Float64,Float64}}`: List of angles (θ, φ) for quadrature.
        - `A_matrix::Matrix{Float64}`: Matrix A used in closure transformation (solves the linear system).
        - `slab_indices::Vector{Int}`: Indices of valid moments in the slab geometry.
        - `target_indices::Vector{Int}`: Indices of moments used for closure.
        - `rotations::AbstractVector{<:AbstractMatrix{Float64}}`: Precomputed rotation matrices for each angle.
    """
    fp_func::F
    Angles::Vector{Tuple{Float64,Float64}}
    A_matrix::Matrix{Float64}
    slab_indices::Vector{Int}
    target_indices::Vector{Int}
    rotations::AbstractVector{<:AbstractMatrix{Float64}}
    function Constants(fp_func::F, Angles::Vector{Tuple{Float64,Float64}}, A_matrix::Matrix{Float64}, slab_indices::Vector{Int}, target_indices::Vector{Int}, rotations::AbstractVector{<:AbstractMatrix{Float64}}) where {F<:Function}
        new{F}(fp_func, Angles, A_matrix, slab_indices, target_indices, rotations)
    end
end

function init_constants(M::Int)
    """
        Initialize constants for Gramian moment equations in 1D with 3D velocity.

        # Arguments
        - `M::Int`: Number of moments.

        # Returns
        - `Constants`: Struct containing precomputed constants.
    """
    # Compile the function
    fp_func = compile_fp(M+1)
    
    # Define Angles - hard-coded for M=4
    angles = [
        (3.14159, 1.5708),
        (0.684719, 4.71239),
        (2.03444, 1.5708),
        (2.18628, 0.886077)
    ]
    
    # Compute Matrix A immediately
    results = [fp_func(theta, phi) for (theta, phi) in angles]
    A_matrix = Matrix(hcat(results...)')
    
    # Get valid indices
    slab_indices = get_valid_indices(M)
    target_indices = nidx(M) # What we use for the closure

    # Precompute rotation matrices
    size_moments_full = length(mainmomindex(M))
    rotations = SVector{4}([SMatrix{size_moments_full, size_moments_full}(rot(M, theta, phi)) for (theta, phi) in angles])

    return Constants(fp_func, angles, A_matrix, slab_indices, target_indices, rotations) # What we use for the closure)
end


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
    constants::Constants

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

        # Setup the Constants
        constants = init_constants(M)

        new{10, n, typeof(Knudsen)}(inv(Knudsen), χ, n, closure_value, constants)
    end
end

# * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = ["U_(000)", "U_(100)", "U_(200)", "U_(020)", "U_(300)", "U_(120)", "U_(400)", "U_(220)", "U_(040)", "U_(022)"]
# ToDo: Change for 3D velocity case
Trixi.varnames(::typeof(cons2prim), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = ntuple(i->"w^($(i-1))", Mp1)

# * basic variable definitions to be used as e.g. sc indicator variable
Trixi.density(u, eqns::GramianMomentEquations1D3D{Mp1}) where {Mp1} = u[1]


"""
    Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
    Flux function: F(U) = u_{k+1}, k=0,…,M where u_{M+1} = C(u0, u1, …, uM)
"""
function Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
    """
        First flux components from shifed moments
        Account for jump skip in indices, when having ∂_t U_200 + ∂_x U_300 = ...
        We have U_200 as index 3, U_300 as index 5, etc.
    """
    known_moments = SVector{6}(ntuple(i->if i <= 2; u[i+1]; else; u[i+2]; end, 6))
    closure_transformation = closure_transform(u, equations)
    return SVector{10}(known_moments..., closure_transformation...)
end


function closure_transform(u, equations)
    """
        Closure transformation for Gramian moment equations in 1D with 3D velocity.

        # Arguments
        - `u`: Input moments in slab geometry.
        - `equations::GramianMomentEquations1D3D`: The equations struct containing constants and closure type.

        # Returns
        - `SVector{4, Float64}`: Transformed moments after applying closure.
    """
    # ToDo: Make generic for arbitrary M
    M = 4

    # Reallocate the moments into the whole geometry with 0 moments for slab
    moments_full = zeros(Float64, length(mainmomindex(M))) # The full set of moments 
    j = 1
    for i in equations.constants.slab_indices
        moments_full[i] = u[j]
        j += 1
    end

    # Initialize rhs and target_indices
    target_indices = equations.constants.target_indices
    rhs = Float64[] # To store the rhs

    # Loop over angles
    for (i, (theta, phi)) in enumerate(equations.constants.Angles)
        # Get Rotation Matrix
        R = equations.constants.rotations[i]
        # Rotate
        rotated_moments_full = R * moments_full

        # Extract (m0, m1, m2, m3, m4)
        substituted = rotated_moments_full[target_indices]

        # Check realizability
        # is_realizable = all(eigen(gramian(substituted, equations.n)).values .>= 0.0)
        # if !is_realizable
        #     println("Warning: Non-realizable moments encountered in closure transformation for angle (θ=$(theta), φ=$(phi)).")
        #     println("Closure transformation for u: ", u)
        #     println("Substituted moments: ", substituted)
        # end
        # Apply closure
        val = closure(substituted, equations)
        push!(rhs, val)
    end

    # Solve for transformed moments
    transformed_moments = equations.constants.A_matrix \ rhs;
    return SVector{4, Float64}(transformed_moments)
end

########################## Closure ##########################

# Overloading to determine Gram/ExtGram/Grad and even/odd case
"""
    closure(u, equations::GramianMomentEquations1D3D)

    Decides which closure implementation is used based on the equations.closure value
"""
closure(u, equations::GramianMomentEquations1D3D) = closure(u, equations, Val(equations.closure))

# -------------------------
# Gramian closure 
# -------------------------
# Decides which closure implementation is used based on the number of Moments (even vs. odd)
# closure(u, equations::GramianMomentEquations1D3D, ::Val{:Gram}) = closure(u, equations, Val(iseven(length(u)-1) ? :GramEven : :GramOdd))
# -------------------------
# closure(u, equations::GramianMomentEquations1D3D, ::Val{:ExtGram}) = closure(u, equations, Val(iseven(length(u)-1) ? :ExtGramEven : :ExtGramOdd))

# Even case (classical and extended)
# ToDo: Nicely rewrite the closures so that they can be used for both 1D and 1D3D cases without code duplication
function closure(u, equations::GramianMomentEquations1D3D, ::Val{:GramEven})
    """
        closure(u, equations::GramianMomentEquations1D3D, ::Val{:GramEven})

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

function closure(u, equations::GramianMomentEquations1D3D, ::Val{:ExtGramEven})
    """
        closure(u, equations::GramianMomentEquations1D3D, ::Val{:ExtGramEven})

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

# todo: implement equilibrium moments for 3D velocity case
# function relaxation_source(u, x, t, equations::GramianMomentEquations1D3D{Mp1}) where {Mp1}
    """
        relaxation_source(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}

        RHS = 1/Kn (u - u_eq)
    """
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
# zero_source(u, x, t, eqns::EqT) where {N, EqT <: Trixi.AbstractEquations{1, N}} = SVector{N}(ntuple(i->0.0, N))


# computes c^n = v^n + 0 + (n over 2) v^(n-2)C^2 + ... + (n over 1) vC^(n-1) + C^n, where u[k] is the primitive moment arising from C^(k-1) 
# function binomial_moment_sum(u, n::Integer=length(u)-1)
#     @assert length(u) > n > 1
#     ρ = u[1]; v = u[2]
#     res = v^n
#     for k=2:n
#         res += binomial(n, k)*v^(n-k)*u[k+1]
#     end
#     return res
# end



# convert primitive [1, v, C^2, ...] variables to conservative [1, c, c^2, ....] variables
# function moment_prim2cons(u_prim)
#     m = length(u_prim)
#     @assert m > 2
#     cons = zeros(eltype(u_prim), m)
#     # u_prim[1] is rho and u_prim[2] is v
#     cons[1] = u_prim[1]; cons[2] = u_prim[1]*u_prim[2]
#     for k=3:m
#         cons[k] = u_prim[1]*binomial_moment_sum(u_prim, k-1)
#     end
#     return cons
# end

# convert conservative to primitive variables
# function moment_cons2prim(u_cons)
#     m = length(u_cons)
#     @assert m > 2
#     prim = zeros(eltype(u_cons), m)
#     prim[1] = u_cons[1]; prim[2] = u_cons[2]/prim[1]
#     for k=3:m
#         prim[k] = u_cons[k]/prim[1]- binomial_moment_sum(prim, k-1)
#     end
#     return prim
# end

# ToDo: Need to check for 3D velocity case
# Trixi.prim2cons(u, eqns::GramianMomentEquations1D3D) = moment_prim2cons(u)
# Trixi.cons2prim(u, eqns::GramianMomentEquations1D3D) = moment_cons2prim(u)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D3D) = u


# Derivative of closure
# function dCdu(u, equations::GramianMomentEquations1D3D)
#     # automatic differentiation
#     closure_wrapped(x) = closure_transform(x, equations)
#     grad = ForwardDiff.gradient(closure_wrapped, u)
#     return ForwardDiff.value(grad)
# end
function dCdu(u, equations::GramianMomentEquations1D3D, h::Float64 = 1e-6)
    """
        Derivative of the closure function w.r.t. the moments u
    """
    grad = zeros(eltype(u), 4, length(u))
    u_aux = zeros(eltype(u), length(u)); @. u_aux = u
    for i in eachindex(u)
        u_aux[i] += h
        grad[:,i] = (closure_transform(u_aux, equations) - closure_transform(u, equations))/h
        u_aux[i] = u[i]
    end
    return grad
end



function flux_jacobian(u, equations::GramianMomentEquations1D3D)
    """
        Jacobian of the flux function
    """
    m = length(u)
    A = zeros(eltype(u), m, m)
    # ToDo: Make generic
    for i=1:2 A[i, i+1] = 1 end
    for i=3:6 A[i, i+2] = 1 end
    A[7:10, :] .= dCdu(u, equations)
    return A
end

function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3D)
    """
        Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
    """
    λ_l = Trixi.max_abs_speeds(u_l, equations)
    λ_r = Trixi.max_abs_speeds(u_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end


function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D3D)
    """
        Estimate the flux Jacobian eigenvalues by means of Gerschgorin
    """
    u = Float64.(u)
    return maximum(abs.(real.(eigen(flux_jacobian(u, equations)).values)))
end




##########################################################
#################### Maximum speeds ######################
##########################################################

# * For now, we assume a constant maximum speed (should be changed later)
# ? What is the optimal maximum speed? We want to avoid calculating the flux Jacobian eigenvalues every time step if possible.
# function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D3D)
#     # REPLACE 1.0 with the maximum microscopic velocity of your model.
#     # For example:
#     # - If your basis is defined on [-1, 1], return 1.0.
#     # - If this is a gas dynamics code, this might need to be higher.
#     # ? What valud to choose here?
#     return 10.0 
# end

# # This remains the same, but now calls the faster constant version above.
# function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3D)
#     λ_l = Trixi.max_abs_speeds(u_l, equations)
#     λ_r = Trixi.max_abs_speeds(u_r, equations)
#     return max(λ_l, λ_r)
# end