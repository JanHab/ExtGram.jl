
struct Maxwellian
    ρ::Real
    v::Real
    θ::Real # defined here as ∫ C^2 f dc (no factor 1/2)
    Maxwellian(ρ, v, θ) = new(ρ, v, θ)
end

(f::Maxwellian)(c::Real) = f.ρ/sqrt(2*π*f.θ) * exp(-(c-f.v)^2 / (2*f.θ)) # note: normalized in 1D velocity space 
# -> ∫f(c)dc = ρ




# Gauss-Hermite quadrature
function convective_moments(f::Maxwellian, ::Val{N}) where {N}
    ξ, w = gausshermite(N+1) # +1 for good measure, should not be necessary
    # compute the integral using the coordinate transform ξ = C/sqrt(2*f.θ)
    # and divide out the exp(-c^2) contained in the gauss hermite quadrature rule (note: quite the inefficient implementation here)
    C = sqrt(2*f.θ) .* ξ; c = C .+ f.v
    fw = f.(c) .* w .* exp.(ξ .^ 2) * sqrt(2*f.θ)
    return SVector{N,Float64}(ntuple(n->sum(c .^(n-1) .* fw), N))
end




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



struct InitialConditionsShockTubeNormalized{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsShockTubeNormalized(f_left, f_right, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        left = convective_moments(f_left, Val(Mp1))
        right = convective_moments(f_right, Val(Mp1))
        # Normalize moments
        # W_n = 1/ρ (1/θ)^{(n)/2} u_n with u_n the n-th convective moment
        # W_n = [1, 0, 1, W_3, …, W_n]'
        left_n = similar(left); right_n = similar(right);
        for i in eachindex(left)
            left_n[i] = 1 / left[1] * (1 / left[3])^((i-1)/2) * left[i]
            right_n[i] = 1 / right[1] * (1 / right[3])^((i-1)/2) * right[i] # same size
        end
        # ToDo: verbose=false
        @assert check_realizability(left_n, verbose=false) && check_realizability(right_n, verbose=false)
        return new{Mp1}(left_n, right_n)
    end
end

function (ic::InitialConditionsShockTubeNormalized)(coords, t, equations::GramianMomentEquations1D)
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end


struct InitialConditionsTwoShocks{N}
    outer::SVector{N}
    inner::SVector{N}

    function InitialConditionsTwoShocks(f_outer, f_inner, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        outer = convective_moments(f_outer, Val(Mp1))
        inner = convective_moments(f_inner, Val(Mp1))
        # ToDo: verbose=false
        @assert check_realizability(outer, verbose=false) && check_realizability(inner, verbose=false)
        return new{Mp1}(outer, inner)
    end
end

function (ic::InitialConditionsTwoShocks)(coords, t, equations::GramianMomentEquations1D)
    if -1.0 < coords[1] && coords[1] < 1.0; return ic.inner; else; return ic.outer; end
end


#######################################
########### 1D3D Maxwellian ###########
#######################################
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

struct InitialConditionsShockTube1D3D{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsShockTube1D3D(f_left, f_right, M, eqns::GramianMomentEquations1D3D{Mp1}) where {Mp1}
        left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
        right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)
        # ToDo: Adapt to 3D velocity case
        # @assert check_realizability(left, verbose=false) && check_realizability(right, verbose=false)
        return new{Mp1}(left, right)
    end
end

function (ic::InitialConditionsShockTube1D3D)(coords, t, equations::GramianMomentEquations1D3D)
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end

#######################################
# Electron hole distribution function #
#######################################
struct ElectronHole
    ψ::Real
    Δ::Real
    v0::Real
    β::Real
    ElectronHole(ψ, Δ, v0, β) = new(ψ, Δ, v0, β)
end

f0_eh(v, ϕ, v0) = 1 / √(2π) * exp(-.5 * (sign(v) * √(v^2 - 2ϕ) - v0)^2)
f1_eh(v, ϕ, v0, β) = 1 / √(2π) * exp(-.5 * (β  * (v^2 - 2ϕ) + v0)^2)
f_eh(v, ϕ, v0, β) = if v^2 > 2ϕ; f0_eh(v, ϕ, v0); else; f1_eh(v, ϕ, v0, β); end;

(eh::ElectronHole)(x::Real, c::Real) = begin
    ϕ = eh.ψ * exp(((x/eh.Δ))^2)
    return f_eh(c, ϕ, eh.v0, eh.β)
end

function convective_moments(eh::ElectronHole, x::Real, ::Val{N}) where {N}
    moments = zeros(N)
    c = LinRange(-5.0, 5.0, 1_000) # todo: maybe make this generic
    # Compute the integral with the trapezoidal rule
    fc = eh.(x, c)
    for n in 1:N
        moments[n] = trapz(c, c.^(n-1) .* fc)
    end
    return SVector{N,Float64}(moments)
end

struct InitialConditionsElectronHole{N}
    eh::ElectronHole
    MP1::Int

    function InitialConditionsElectronHole(eh, eqns::GramianMomentEquations1D{Mp1}) where {Mp1}
        return new{Mp1}(eh, Mp1)
    end
end

function (ic::InitialConditionsElectronHole)(coords, t, equations::GramianMomentEquations1D)
    moments = convective_moments(ic.eh, coords[1], Val(ic.MP1))
    # ToDo: verbose=false
    @assert check_realizability(moments, verbose=false)
    return moments
end