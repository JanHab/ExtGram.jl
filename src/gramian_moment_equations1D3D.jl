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

function init_constants(M::Int, angles::Vector{Tuple{Float64,Float64}} = [
        (3.14159, 1.5708),
        (0.684719, 4.71239),
        (2.03444, 1.5708),
        (2.18628, 0.886077)
    ])
    """
        Initialize constants for Gramian moment equations in 1D with 3D velocity.

        # Arguments
        - `M::Int`: Number of moments.
        - `angles::Vector{Tuple{Float64,Float64}}`: List of angles (θ, φ) for quadrature.

        # Returns
        - `Constants`: Struct containing precomputed constants.
    """
    # Compile the function
    fp_func = compile_fp(M+1)
    
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
struct GramianMomentEquations1D3D{Mp1, N, RealT <: Real} <: GramianMomentEquations{Mp1, N, RealT}
    inv_Kn::RealT # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    n::Int
    closure::Symbol
    constants::Constants

    function GramianMomentEquations1D3D(M::Integer, Knudsen::Real, closure::String="ExtGram"; χ_set="optimal", angles::Vector{Tuple{Float64,Float64}} = [
        (3.14159, 1.5708),
        (0.684719, 4.71239),
        (2.03444, 1.5708),
        (2.18628, 0.886077)
    ])
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
        elseif closure == "Grad"
            closure_value = :Grad
        else
            error("Unknown closure type: $closure. Supported types are \"Gram\" and \"ExtGram\".\n.")
        end

        # Setup the Constants
        constants = init_constants(M, angles)

        new{10, n, typeof(Knudsen)}(inv(Knudsen), χ, n, closure_value, constants)
    end
end

# * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = (
    "U_(000)", "U_(100)", "U_(200)", "U_(020)", "U_(300)", 
    "U_(120)", "U_(400)", "U_(220)", "U_(040)", "U_(022)"
)
Trixi.varnames(::typeof(cons2prim), ::GramianMomentEquations1D3D{Mp1}) where {Mp1} = (
    "W_(000)", "W_(100)", "W_(200)", "W_(020)", "W_(300)", 
    "W_(120)", "W_(400)", "W_(220)", "W_(040)", "W_(022)"
)

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

        # Apply closure
        val = closure(substituted, equations)
        push!(rhs, val)
    end

    # Solve for transformed moments
    transformed_moments = equations.constants.A_matrix \ rhs;
    return SVector{4, Float64}(transformed_moments)
end


##########################################################
##################### Source term ########################
##########################################################
function zero_source(u, x, t, equations::GramianMomentEquations1D3D)
    return SVector{10}(ntuple(i->0.0, 10))
end
function relaxation_source(u, x, t, equations::GramianMomentEquations1D3D)
    """
        relaxation_source(u, x, t, equations::GramianMomentEquations1D3D)

        RHS = -1/Kn (u - u_eq)
    """
    # 1. Calculate Physical Properties
    rho = u[1]
    
    # Calculate Central Moments P (needed for Temperature)
    p_prim = cons2prim(u, equations)
    
    P_200 = p_prim[3] # P_xx
    P_020 = p_prim[4] # P_yy
    
    # Slab Geometry: Transverse symmetry P_zz = P_yy (P_002 = P_020)
    P_002 = P_020 
    
    # Temperature theta = (U_200 + U_020 + U_002) / (3 * rho)
    theta = (P_200 + P_020 + P_002) / (3 * rho)
    
    # 2. Compute Equilibrium Central Moments P_eq
    # P_eq_{ijk} = rho * Gauss(i) * Gauss(j) * Gauss(k)
    p_eq = SVector{10, Float64}(ntuple(idx -> begin
        i, j, k = MOMENT_INDICES_1D3D[idx] # indices for (i, j, k)
        
        # If any power is odd, central moment equilibrium is 0
        if isodd(i) || isodd(j) || isodd(k)
            return 0.0
        end
        
        # Gaussian moments: (n-1)!! * theta^(n/2)
        val_i = double_factorial(i-1) * theta^(i/2)
        val_j = double_factorial(j-1) * theta^(j/2)
        val_k = double_factorial(k-1) * theta^(k/2)
        
        return rho * val_i * val_j * val_k
    end, 10))
    
    # 3. Transform P_eq back to Conservative U_eq
    u_eq = prim2cons(p_eq, equations)
    
    # 4. Return DVM source
    return SVector{10}(-equations.inv_Kn .* (u .- u_eq))
end
# Not necessary to redefine zero_source, as the default implementation for the 1D1D case suffices.

# 1. Define the moment structure for the 10-moment system
# Mapping indices 1..10 to (i, j, k) powers
const MOMENT_INDICES_1D3D = (
    (0,0,0), # 1: rho
    (1,0,0), # 2: rho * v
    (2,0,0), # 3: P_xx
    (0,2,0), # 4: P_yy
    (3,0,0), # 5
    (1,2,0), # 6
    (4,0,0), # 7
    (2,2,0), # 8
    (0,4,0), # 9
    (0,2,2)  # 10: P_yyzz (mixed transverse)
)

function double_factorial(n::Int)
    """
        # Helper for double factorial (n-1)!!

        (n-1)!! = (n-1) * (n-3) * (n-5) * ... until 1 or 2
    """
    n <= 0 && return 1.0
    val = 1.0
    for k in 1:2:n
        val *= k
    end
    return val
end
# double_factorial(n) = prod((n):-2:1)

@inline function get_moment_val(u, i_req, j_req, k_req, indices_list)
    """
        get_moment_val(u, i, j, k, indices_list)

    Helper to retrieve U_{ijk} from the state vector `u` based on the provided 
    indices_list. Returns 0.0 if the moment is not in the system.
    """
    # Search for the tuple in the index list (constant folding should optimize this)
    for (idx, (i, j, k)) in enumerate(indices_list)
        if i == i_req && j == j_req && k == k_req
            return u[idx]
        end
    end
    return 0.0 # Return 0 if moment not tracked (or handle error)
end


function moment_cons2prim(u_cons, eqns::GramianMomentEquations1D3D, moment_indices=MOMENT_INDICES_1D3D)
    # Type and length of array
    T = eltype(u_cons)
    N = length(moment_indices)

    # Initialize Mutable Vector
    prim = MVector{N, T}(undef)

    # Explicitly set the first two moments (Avoiding the if-check)
    prim[1] = u_cons[1]; prim[2] = u_cons[2] / prim[1] # rho, v

    # Loop for the remaining moments (Indices 3 to N)
    v = prim[2]
    for idx in 3:N
        i, j, k = moment_indices[idx]
        val = zero(T)
        
        # Calculate P_{ijk} = sum Binomial * (-v)^... * U
        for m in 0:i
            u_mjk = get_moment_val(u_cons, m, j, k, moment_indices)
            val += binomial(i, m) * (-v)^(i-m) * u_mjk
        end
        
        prim[idx] = val
    end

    return SVector(prim)
end


function moment_prim2cons(u_prim, eqns::GramianMomentEquations1D3D, moment_indices=MOMENT_INDICES_1D3D)
    T = eltype(u_prim)
    N = length(moment_indices)

    # Use MVector
    cons = MVector{N, T}(undef)

    # Explicitly set the first two conservative moments
    cons[1] = u_prim[1]; cons[2] = u_prim[1]*u_prim[2]  # rho, rho*v

    # Loop for the rest
    v = u_prim[2]
    for idx in 3:N
        i, j, k = moment_indices[idx]
        val = zero(T)
        for m in 0:i
            # Inner check: ensure we use 0.0 for P_{100}, not 'v'
            # ? Is this necessary? Do we even get here? Is it correct?
            if m == 1 && j == 0 && k == 0
                p_mjk = 0.0
            else
                p_mjk = get_moment_val(u_prim, m, j, k, moment_indices)
            end
            
            val += binomial(i, m) * v^(i-m) * p_mjk
        end
        cons[idx] = val
    end

    return SVector(cons)
end

# Link to Trixi
Trixi.cons2prim(u, eqns::GramianMomentEquations1D3D) = moment_cons2prim(u, eqns)
Trixi.prim2cons(u, eqns::GramianMomentEquations1D3D) = moment_prim2cons(u, eqns)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D3D) = u


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




##########################################################
#################### Maximum speeds ######################
##########################################################

# This remains the same, but now calls the faster constant version above.
function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3D)
    """
        λ_max = v_x ± γ √(θ), with γ=5
    """
    ρ_l = u_l[1]; ρ_r = u_r[1]

    # Calculate Central Moments P (needed for Temperature)
    p_prim_l = cons2prim(u_l, equations); p_prim_r = cons2prim(u_r, equations);
    
    P_200_l = p_prim_l[3]; P_200_r = p_prim_r[3] # P_xx
    P_020_l = p_prim_l[4]; P_020_r = p_prim_r[4] # P_yy
    
    # Slab Geometry: Transverse symmetry P_zz = P_yy (P_002 = P_020)
    P_002_l = P_020_l; P_002_r = P_020_r 

    v_kl = p_prim_l[2]; v_kr = p_prim_r[2] # velocity in x-direction (slab geometry)
    
    # Temperature theta = (U_200 + U_020 + U_002) / (3 * rho)
    θ_l = 1 / (3 * ρ_l) * (P_200_l + P_020_l + P_002_l - ρ_l * v_kl^2)
    θ_r = 1 / (3 * ρ_r) * (P_200_r + P_020_r + P_002_r - ρ_r * v_kr^2)
    
    γ = 5.0
    λ_l = ρ_l + γ * sqrt(θ_l)
    λ_r = ρ_r + γ * sqrt(θ_r)
    return max(λ_l, λ_r)
end