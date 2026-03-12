
"""
    Maxwellian distribution function in 1D velocity space
"""
struct Maxwellian
    ρ::Real
    v::Real
    θ::Real # defined here as ∫ C^2 f dc (no factor 1/2)
    Maxwellian(ρ, v, θ) = new(ρ, v, θ)
end

(f::Maxwellian)(c::Real) = f.ρ/sqrt(2*π*f.θ) * exp(-(c-f.v)^2 / (2*f.θ)) # note: normalized in 1D velocity space 
# -> ∫f(c)dc = ρ




# Gauss-Hermite quadrature
"""
    Compute the integral using the coordinate transform ξ = C/sqrt(2*f.θ) and divide out the exp(-c^2) contained in the gauss hermite quadrature rule (note: quite the inefficient implementation here)
"""
function convective_moments(f::Maxwellian, ::Val{N}) where {N}
    ξ, w = gausshermite(N+1) # +1 for good measure, should not be necessary
    C = sqrt(2*f.θ) .* ξ; c = C .+ f.v
    fw = f.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*f.θ)
    return SVector{N,Float64}(ntuple(n->sum(c .^(n-1) .* fw), N))
end



"""
    Compute the convective moments ∫ C^{n-1} f dc where C = c - v
    and return the primitive moments [ρ, v, C^2, C^3, ..., C^{N-1}]
"""
function primitive_moments(f::Maxwellian, ::Val{N}) where {N}
    ξ, w = gausshermite(N+1)
    C = sqrt(2*f.θ) .* ξ; c = C .+ f.v
    fw = f.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*f.θ)
    moments = ntuple(n->sum(C .^(n-1) .* fw), N)
    if N > 1
        @assert abs.(moments[2]) < 1e-12
        # divide density out of higher moments (e.g. C^2 gives ρθ -> remove ρ)
        return SVector{N, Float64}(moments[1], f.v, (moments[3:end] ./ moments[1])...)
    else
        return SVector{N, Float64}(moments[1])
    end
end



"""
    Shock tube initial conditions with left and right distribution functions
"""
struct InitialConditionsShockTube{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsShockTube(f_left, f_right, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        left = convective_moments(f_left, Val(Mp1))
        right = convective_moments(f_right, Val(Mp1))
        # ToDo: verbose=false
        @assert check_realizability(left, verbose=false) && check_realizability(right, verbose=false)
        return new{Mp1}(left, right)
    end
end

function (ic::InitialConditionsShockTube)(coords, t, equations::GramianMomentEquations1D)
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end




#######################################
########### 1D3D Maxwellian ###########
#######################################
"""
    Maxwellian distribution function in 3D velocity space (1D spatial)
"""
struct Maxwellian1D3D
    ρ::Real
    v::NTuple{3,Real}
    θ::Real # defined here as ∫ C^2 f dc (no factor 1/2)
    Maxwellian1D3D(ρ, v, θ) = new(ρ, v, θ)
end

(f::Maxwellian1D3D)(cx::Real, cy::Real, cz::Real) = f.ρ/(2*π*f.θ)^(3/2) * exp(-((cx-f.v[1])^2 + (cy-f.v[2])^2 + (cz-f.v[3])^2) / (2*f.θ)) # note: normalized in 1D velocity space 

"""
    index3D(total_degree)

    Return a vector of multi-indices (i,j,k) with i+j+k == total_degree.
    Ordering matches a lexicographic-style ordering and is intended for moment lists.
"""
function index3D(n::Integer)
    tuples = collect(Iterators.product(0:n, 0:n, 0:n))
    vecs = [collect(t) for t in tuples]
    selected = filter(v -> sum(v) == n, vecs)
    sort!(selected)
    return reverse(selected)
end

"""
    multi_index_list(M)

    Return a flattened list of multi-indices for all shells 0..M.
"""
function multi_index_list(M::Integer)
    return vcat([index3D(n) for n in 0:M]...)
end

"""
    convective_moments_1D3D(M; rho=1.0, u=(0,0,0), theta=1.0, q=8)

    Compute the velocity moments up to `M` (all multi-indices with i+j+k <= M)
    using tensor-product Gauss–Hermite with `q` nodes per dimension. Returns a Vector{Float64}
    with the same ordering as `multi_index_list(M)`.
"""
function convective_moments_1D3D(M::Integer, rho::Real=1.0, u::NTuple{3,Real}=(0.0,0.0,0.0), theta::Real=1.0; q::Integer=35+1)
    # 35+1 is max degree in the test cases for M=4 (hard-coded)
    # quadrature nodes and weights for ∫ e^{-x^2} g(x) dx
    ξ, w = gausshermite(q)

    # scaling: v = u + sqrt(2 theta) * ξ
    s = sqrt(2 * theta)
    ux, uy, uz = u

    # Precompute 1D transformed node values for each dimension
    vx_nodes = ux .+ s .* ξ
    vy_nodes = uy .+ s .* ξ
    vz_nodes = uz .+ s .* ξ

    # Prefactor from change of variables: rho / π^{3/2}
    pref = rho / (π^(3/2))

    # Prepare storage for moments: a Dict from (i,j,k) -> value
    idxs = multi_index_list(M)
    moments = zeros(Float64, length(idxs))

    # Triple loop over quadrature points
    # Complexity q^3 but q is typically small (8..16)
    for a in 1:q, b in 1:q, c in 1:q
        wt = w[a] * w[b] * w[c]
        vx = vx_nodes[a]
        vy = vy_nodes[b]
        vz = vz_nodes[c]
        # Accumulate into all requested moments
        for (idx_pos, (i,j,k)) in enumerate(idxs)
            moments[idx_pos] += wt * (vx^i) * (vy^j) * (vz^k)
        end
    end

    # Apply prefactor
    moments .*= pref
    return moments
end

"""
    Shock tube initial conditions with left and right distribution functions in 1D3D
"""
struct InitialConditionsShockTube1D3D{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsShockTube1D3D(f_left, f_right, M, eqns::GramianMomentEquations1D3D{Mp1}) where {Mp1}

        # Compute full convective moments
        left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
        right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)

        # Reduce to slab indices only
        slab_indices = get_valid_indices(M)
        left_reduced = SVector{length(slab_indices)}(left[slab_indices])
        right_reduced = SVector{length(slab_indices)}(right[slab_indices])

        return new{length(slab_indices)}(left_reduced, right_reduced)
    end
end

function (ic::InitialConditionsShockTube1D3D)(coords, t, equations::GramianMomentEquations1D3D)
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end




#############################################
########## Vlasov Example Initial Conditions ##########
#############################################

"""
        Initial conditions for Landau damping problem in Vlasov-Poisson system

        The probability distribution function is a small perturbtation of a Maxwellian:
        f(x,c) = 1 / √(2π) * (1 + ϵ * cos(k * x)) * exp(- c^2 / 2)
    """
struct InitialConditionsLandauDamping{N}
    ρ0::Float64
    ϵ::Float64
    v0::Float64
    θ0::Float64
    k::Float64

    function InitialConditionsLandauDamping(ρ0::Float64, ϵ::Float64, v0::Float64, θ0::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(ρ0, ϵ, v0, θ0, k)
    end
end

function (ic::InitialConditionsLandauDamping)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    ρx = ic.ρ0 * (1 + ic.ϵ * cos(ic.k * coords[1]))
    f = Maxwellian(ρx, ic.v0, ic.θ0)
    return convective_moments(f, Val(Mp1))
end

"""
    Initial conditions for the two-stream instability problem in Vlasov-Poisson system

    The probability distribution function is given by:
        f(x,c) = 1 / √(2π) * (1 + ϵ cos(k * x)) * (exp(- c^2 / 2) * c^2
"""
struct InitialConditionsTwoStream{N}
    ϵ::Float64
    k::Float64

    function InitialConditionsTwoStream(ϵ::Float64, k::Float64, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(ϵ, k)
    end
end

function (ic::InitialConditionsTwoStream)(coords, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    max(c) = 1/sqrt(2*π) * exp(-c^2 / 2) * c^2 .* (1 + ic.ϵ * cos(ic.k * coords[1])) # note: normalized in 1D velocity space 
    # copying from convective_moments for standard Maxwellian
    ξ, w = gausshermite(Mp1+1) # +1 for good measure, should not be necessary
    C = sqrt(2*1.0) .* ξ; c = C .+ 0.0
    fw = max.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*1.0)
    return SVector{Mp1,Float64}(ntuple(n->sum(c .^(n-1) .* fw), Mp1))
end