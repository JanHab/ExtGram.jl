
struct Maxwellian
    """
        Maxwellian distribution function in 1D velocity space
    """
    ρ::Real
    v::Real
    θ::Real # defined here as ∫ C^2 f dc (no factor 1/2)
    Maxwellian(ρ, v, θ) = new(ρ, v, θ)
end

(f::Maxwellian)(c::Real) = f.ρ/sqrt(2*π*f.θ) * exp(-(c-f.v)^2 / (2*f.θ)) # note: normalized in 1D velocity space 
# -> ∫f(c)dc = ρ




# Gauss-Hermite quadrature
function convective_moments(f::Maxwellian, ::Val{N}) where {N}
    """
        Compute the integral using the coordinate transform ξ = C/sqrt(2*f.θ) and divide out the exp(-c^2) contained in the gauss hermite quadrature rule (note: quite the inefficient implementation here)
    """
    ξ, w = gausshermite(N+1) # +1 for good measure, should not be necessary
    C = sqrt(2*f.θ) .* ξ; c = C .+ f.v
    fw = f.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*f.θ)
    return SVector{N,Float64}(ntuple(n->sum(c .^(n-1) .* fw), N))
end




function primitive_moments(f::Maxwellian, ::Val{N}) where {N}
    """
        Compute the convective moments ∫ C^{n-1} f dc where C = c - v
        and return the primitive moments [ρ, v, C^2, C^3, ..., C^{N-1}]
    """
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




struct InitialConditionsShockTube{N}
    """
        Shock tube initial conditions with left and right distribution functions
    """
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
struct Maxwellian1D3D
    """
        Maxwellian distribution function in 3D velocity space (1D spatial)
    """
    ρ::Real
    v::NTuple{3,Real}
    θ::Real # defined here as ∫ C^2 f dc (no factor 1/2)
    Maxwellian1D3D(ρ, v, θ) = new(ρ, v, θ)
end

(f::Maxwellian1D3D)(cx::Real, cy::Real, cz::Real) = f.ρ/(2*π*f.θ)^(3/2) * exp(-((cx-f.v[1])^2 + (cy-f.v[2])^2 + (cz-f.v[3])^2) / (2*f.θ)) # note: normalized in 1D velocity space 

function index3D(n::Integer)
    """
        index3D(total_degree)

    Return a vector of multi-indices (i,j,k) with i+j+k == total_degree.
    Ordering matches a lexicographic-style ordering and is intended for moment lists.
    """
    tuples = collect(Iterators.product(0:n, 0:n, 0:n))
    vecs = [collect(t) for t in tuples]
    selected = filter(v -> sum(v) == n, vecs)
    sort!(selected)
    return reverse(selected)
end

function multi_index_list(M::Integer)
    """
        multi_index_list(M)

    Return a flattened list of multi-indices for all shells 0..M.
    """
    return vcat([index3D(n) for n in 0:M]...)
end

function convective_moments_1D3D(M::Integer, rho::Real=1.0, u::NTuple{3,Real}=(0.0,0.0,0.0), theta::Real=1.0; q::Integer=35+1)
    """
        convective_moments_1D3D(M; rho=1.0, u=(0,0,0), theta=1.0, q=8)

    Compute the velocity moments up to `M` (all multi-indices with i+j+k <= M)
    using tensor-product Gauss–Hermite with `q` nodes per dimension. Returns a Vector{Float64}
    with the same ordering as `multi_index_list(M)`.
    """
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

struct InitialConditionsShockTube1D3D{N}
    """
        Shock tube initial conditions with left and right distribution functions in 1D3D
    """
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
