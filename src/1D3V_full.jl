"""
    Implementation of the one-dimensional Gramian moment equations with full 3D velocity
    dependence ("1D-3V"), i.e. **without** assuming slab symmetry in the transverse
    directions.

    The system evolves every moment U_(αβγ) = ∫ f v_x^α v_y^β v_z^γ dv with α+β+γ ≤ M,

        ∂_t U_(αβγ) + ∂_x U_(α+1,βγ) = S_(αβγ),

    so the flux of an equation is again a moment of the system, except for the moments of
    total order M+1 that leave it. Those are supplied by the closure.

    # Closure

    The (extended) Gramian closure is a genuinely one-dimensional construction. It is
    lifted to 3D velocity space by evaluating it along a set of directions n(θ,φ):
    rotating the state into a frame whose x-axis is n turns the moments U_(k00) of the
    rotated state into the directional moments ∫ f (v·n)^k dv, k = 0…M, which the 1D
    closure maps to ∫ f (v·n)^(M+1) dv. Expanding that contraction,

        ∫ f (v·n)^(M+1) dv = Σ_{α+β+γ=M+1} c_(αβγ)(n) U_(αβγ),

    gives one linear equation per direction for the order-(M+1) moments. Note that the sum
    runs over **all** order-(M+1) multi-indices, including the α = 0 ones, even though only
    the α ≥ 1 moments appear in the x-flux. The α = 0 columns carry O(1) coefficients, so
    they cannot be dropped: the full shell of

        N_shell = (M+2)(M+3)/2

    moments is solved for, using N_shell directions, and the α ≥ 1 prefix is handed to the
    flux. (The multi-index ordering of `index(M+1)` is descending, so the α ≥ 1 moments are
    exactly the first N_closure = (M+1)(M+2)/2 entries -- no permutation is needed.)

    # Choice of directions

    A rank-(M+1) contraction with M+1 odd satisfies c(-n) = -c(n), so antipodal directions
    carry no independent information and any set containing an antipodal pair is singular.
    The default set is therefore a Fibonacci spiral on the **hemisphere**, which is far
    better conditioned than random directions (cond ≈ 16 vs. ≈ 7·10³ median for M = 4).
"""

"""
    U_t_index(nmax, slab_geometry=false)

    Multi-indices of the evolved moments, i.e. the left-hand sides ∂_t U_(αβγ).
"""
function U_t_index(nmax::Integer, slab_geometry::Bool=false)
    if slab_geometry
        # ToDo: Rename get_valid_indices to something more descriptive, e.g., get_slab_indices
        return mainmomindex(nmax)[get_valid_indices(nmax)]
    else
        return mainmomindex(nmax)
    end
end

"""
    U_x_index(nmax, slab_geometry=false)

    Multi-indices of the flux moments, i.e. U_(αβγ) shifted by one in x.
"""
function U_x_index(nmax::Integer, slab_geometry::Bool=false)
    U_t = U_t_index(nmax, slab_geometry)
    return [U_t[i] .+ [1, 0, 0] for i in 1:length(U_t)]
end

"""
    slab_scatter_matrix(full_indices, slab_indices)

    A moment with odd β or γ vanishes (zero row), every other moment equals its representative
    U_(α, max(β,γ), min(β,γ)). The returned matrix S is the linear map from the full moment vector to the slab-symmetric one, i.e. S * u_slab = u_full
"""
function slab_scatter_matrix(full_indices::Vector{Vector{Int}}, slab_indices::Vector{Vector{Int}})
    position = Dict(ix => j for (j, ix) in enumerate(slab_indices))
    S = zeros(Float64, length(full_indices), length(slab_indices))
    for (p, ix) in enumerate(full_indices)
        (iseven(ix[2]) && iseven(ix[3])) || continue
        S[p, position[[ix[1], max(ix[2], ix[3]), min(ix[2], ix[3])]]] = 1.0
    end
    return S
end


##########################################################
################ Direction (angle) helpers ###############
##########################################################

"""
    direction_from_angles(theta, phi)

    The unit vector n along which the 1D closure is evaluated for the angles (θ, φ).

    This is the first row of `get_rotation_matrix(theta, phi)`, which is the only row the
    closure ever reads: the (k,0,0) component of a rank-k tensor transformation is the
    product of k entries of that row.
"""
direction_from_angles(theta::Real, phi::Real) =
    SVector(-cos(theta) * sin(phi), -sin(theta) * sin(phi), cos(phi))

"""
    angles_from_direction(n)

    Inverse of [`direction_from_angles`](@ref): the (θ, φ) reproducing the unit vector `n`.
"""
function angles_from_direction(n)
    phi = acos(clamp(n[3], -one(eltype(n)), one(eltype(n))))
    s = sin(phi)
    abs(s) < 1e-12 && return (zero(phi), phi)   # n ∥ e_z, θ is arbitrary
    return (atan(-n[2] / s, -n[1] / s), phi)
end

"""
    fibonacci_hemisphere_angles(n_points)

    `n_points` well-separated directions on the upper hemisphere, as (θ, φ) vectors.

    The hemisphere (rather than the full sphere) matters: for odd tensor rank the
    contraction is odd in n, so antipodal directions duplicate each other and a full-sphere
    spiral produces a numerically singular system.
"""
function fibonacci_hemisphere_angles(n_points::Integer)
    golden_angle = π * (3 - sqrt(5.0))
    theta = Vector{Float64}(undef, n_points)
    phi = Vector{Float64}(undef, n_points)
    for i in 0:(n_points - 1)
        z = (i + 0.5) / n_points
        r = sqrt(max(0.0, 1 - z^2))
        a = golden_angle * i
        theta[i + 1], phi[i + 1] = angles_from_direction((r * cos(a), r * sin(a), z))
    end
    return theta, phi
end


##########################################################
################# Equations definition ###################
##########################################################

"""
    GramianMomentEquations1D3V{Mp1, N, MC, NC, NS, RealT}

    Gramian moment equations in one space dimension with full 3D velocity dependence.

    # Type parameters
    - `Mp1`: number of equations (Trixi's NVARS) -- despite the inherited name this is
      `binomial(M+3, 3)`, **not** `M+1`. The abstract supertype fixes the position.
    - `N`  : closure order `n` (see [`closure`](@ref)).
    - `MC` : `M+1`, the number of directional moments the 1D closure consumes.
    - `NC` : number of closing moments the flux needs, `(M+1)(M+2)/2`.
    - `NS` : size of the order-(M+1) shell that is solved for, `(M+2)(M+3)/2`.

    # Constructor
        GramianMomentEquations1D3V(
        M, Knudsen, closure="ExtGram";
        theta=Float64[], phi=Float64[]
    )

    `theta`/`phi` must hold `NS` directions; if omitted, a Fibonacci hemisphere is used.
"""
struct GramianMomentEquations1D3V{Mp1, N, MC, NC, NS, RealT <: Real} <: GramianMomentEquations{Mp1, N, RealT}
    inv_Kn::RealT   # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    M::Int          # highest moment order
    n::Int
    N_equations::Int
    N_closure::Int  # number of moments the closure has to supply to the flux
    N_shell::Int    # number of order-(M+1) moments solved for (N_closure of which are used)
    closure::Symbol
    slab_geometry::Bool
    _U_t_index::Vector{Vector{Int}}
    _U_x_index::Vector{Vector{Int}}
    # Flux lookup table, one entry per equation i:
    #   j > 0  ->  U_(_U_x_index[i]) is itself an evolved moment, found at u[j]
    #   j < 0  ->  U_(_U_x_index[i]) is a closing moment, found at closure_moments[-j]
    _flux_index::NTuple{Mp1, Int}
    # Multi-indices of the closing moments, in the order the closure returns them
    _closure_index::Vector{Vector{Int}}
    # Per direction: the (k,0,0) rows of the rank-≤M rotation, so that _P[i] * u are the
    # directional moments ∫ f (v·n_i)^k dv, k = 0…M, that the 1D closure consumes.
    _P::Vector{Matrix{RealT}}
    # Rows 1:N_closure of the inverse order-(M+1) shell transformation (NC × NS)
    _T_closure::Matrix{RealT}
    # Directions of the closure evaluations
    theta::Vector{RealT}
    phi::Vector{RealT}
    # (α,β,γ) powers of the evolved moments, as tuples (cheap to destructure)
    _pow::Vector{NTuple{3, Int}}
    # _pos[α+1, β+1, γ+1] = position of U_(αβγ) in the state, or 0 if not evolved
    _pos::Array{Int, 3}
    # Positions of P_(200), P_(020), P_(002) in the primitive vector (for the temperature)
    _pressure_index::NTuple{3, Int}

    function GramianMomentEquations1D3V(
        M::Integer, Knudsen::Real, closure::String="ExtGram";
        slab_geometry::Bool=false,
        theta::Vector{<:Real}=Float64[], phi::Vector{<:Real}=Float64[],
    )
        @assert M > 1 "M must be greater than 1, got M = $M."
        @assert length(theta) == length(phi) "theta and phi must have equal length, got $(length(theta)) and $(length(phi))."

        RealT = float(typeof(Knudsen))

        if iseven(M)
            n = Int(M / 2)
            χ = (n + 1) / n
        else
            n = Int((M + 1) / 2)
            χ = (n + 1) / (2n)
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

        _U_t_index = U_t_index(M, slab_geometry)
        _U_x_index = U_x_index(M, slab_geometry)
        N_equations = length(_U_t_index)

        # Resolve every flux moment U_(α+1, β, γ) against the evolved moments U_(α, β, γ).
        # Everything of total order ≤ M is found in _U_t_index; the shifted moments of the
        # order-M block leave the system (total order M+1) and have to come from the closure.
        position_in_u = Dict(index => i for (i, index) in enumerate(_U_t_index))
        flux_index = Vector{Int}(undef, N_equations)
        _closure_index = Vector{Vector{Int}}()
        for i in 1:N_equations
            j = get(position_in_u, _U_x_index[i], 0)
            if j != 0
                flux_index[i] = j
            else
                push!(_closure_index, _U_x_index[i])
                flux_index[i] = -length(_closure_index)
            end
        end
        N_closure = length(_closure_index)

        # The order-(M+1) shell. _closure_index is exactly its α ≥ 1 prefix, because
        # index(M+1) is ordered by descending multi-index.
        shell = slab_geometry ? index(M + 1)[index_1d(M + 1)] : index(M + 1) # ?index(M + 1)
        N_shell = length(shell)
        @assert _closure_index == shell[1:N_closure] "closure moments are not the leading block of index(M+1); the prefix assumption in closure_moments() is violated."
        # ? shell = index(M)
        # ? N_shell = length(shell)

        # Directions: default to a Fibonacci hemisphere (see the module docstring above).
        if isempty(theta)
            theta, phi = fibonacci_hemisphere_angles(N_shell)
        end
        @assert length(theta) == N_shell "the closure needs exactly N_shell = $N_shell directions for M = $M (the full order-$(M+1) shell), got $(length(theta))."

        N_full = length(mainmomindex(M))
        S_state = slab_geometry ? slab_scatter_matrix(mainmomindex(M), _U_t_index) : Matrix{Float64}(I, N_full, N_full)
        S_shell = slab_geometry ? slab_scatter_matrix(index(M + 1), shell) : Matrix{Float64}(I, N_shell, N_shell)

        # Per-direction rotation of the evolved moments onto the directional moments
        # ∫ f (v·n)^k dv, k = 0…M, which are located at rows nidx(M) of the full rotation.
        target_rows = nidx(M)
        _P = [Matrix{RealT}(rot(M, theta[i], phi[i])[target_rows, :] * S_state) for i in 1:N_shell]

        # Linear map from the order-(M+1) shell onto the directional closure values
        T_shell = Matrix{RealT}(undef, N_shell, N_shell)
        for i in 1:N_shell
            T_shell[i, :] = S_shell' * tensor_transformation(M + 1, theta[i], phi[i])[1, :]
            # ? T_shell[i, :] = tensor_transformation(M, theta[i], phi[i])[1, :]
        end
        κ = cond(T_shell)
        if κ > 1e8
            @warn "The closure direction set is badly conditioned (cond = $κ); the recovered order-$(M+1) moments will be inaccurate. Antipodal directions are degenerate for odd M+1 -- keep all directions on one hemisphere." M κ
        end
        # Only the α ≥ 1 prefix of the solution is ever used, so keep only those rows.
        _T_closure = Matrix{RealT}(inv(T_shell)[1:N_closure, :])

        # Fast (α,β,γ) -> position lookup, and the powers as tuples
        _pow = [(ix[1], ix[2], ix[3]) for ix in _U_t_index]
        _pos = zeros(Int, M + 1, M + 1, M + 1)
        for (p, ix) in enumerate(_U_t_index)
            _pos[ix[1] + 1, ix[2] + 1, ix[3] + 1] = p
        end
        if slab_geometry
            # Mirror-aware lookup: a non-evolved moment with β, γ even equals its
            # representative U_(α, max(β,γ), min(β,γ)); moments with an odd power vanish
            # by symmetry and keep position 0. moment_val / moment_or_zero / temperature
            # (via _pressure_index below) then work in slab mode without any changes.
            for ix in mainmomindex(M)
                (iseven(ix[2]) && iseven(ix[3])) || continue
                _pos[ix[1] + 1, ix[2] + 1, ix[3] + 1] =
                    _pos[ix[1] + 1, max(ix[2], ix[3]) + 1, min(ix[2], ix[3]) + 1]
            end
        end
        # cons2prim/prim2cons rely on ρ and ρ v_x occupying the first two slots (in full
        # mode additionally ρ v_y, ρ v_z in slots 3:4).
        if slab_geometry
            @assert _pow[1] == (0, 0, 0) && _pow[2] == (1, 0, 0) "unexpected moment ordering; slab mode assumes U_(000), U_(100) in slots 1:2."
        else
            @assert _pow[1] == (0, 0, 0) && _pow[2] == (1, 0, 0) && _pow[3] == (0, 1, 0) && _pow[4] == (0, 0, 1) "unexpected moment ordering; cons2prim/prim2cons assume U_(000), U_(100), U_(010), U_(001) in slots 1:4."
        end
        _pressure_index = (_pos[3, 1, 1], _pos[1, 3, 1], _pos[1, 1, 3])

        # NOTE: the first type parameter is Trixi's NVARS (see the abstract type in
        # src/gramian_moment_equations.jl), so it must be N_equations, not M+1.
        new{N_equations, n, M + 1, N_closure, N_shell, RealT}(
            RealT(inv(Knudsen)), RealT(χ), Int(M), n, N_equations, N_closure, N_shell,
            closure_value, slab_geometry, _U_t_index, _U_x_index, Tuple(flux_index), _closure_index,
            _P, _T_closure,
            convert(Vector{RealT}, theta), convert(Vector{RealT}, phi),
            _pow, _pos, _pressure_index)
    end
end

# * Conservative to Primitive variables and vice versa
_moment_name(prefix, ix) = string(prefix, "_(", ix[1], ix[2], ix[3], ")")

Trixi.varnames(::typeof(cons2cons), equations::GramianMomentEquations1D3V{Mp1}) where {Mp1} =
    ntuple(i -> _moment_name("U", equations._U_t_index[i]), Val(Mp1))
# The primitive variables are the *central* moments of the same multi-indices (with the
# three velocity components in slots 2:4), not the flux moments.
Trixi.varnames(::typeof(cons2prim), equations::GramianMomentEquations1D3V{Mp1}) where {Mp1} =
    ntuple(i -> _moment_name("W", equations._U_t_index[i]), Val(Mp1))

# * basic variable definitions to be used as e.g. sc indicator variable
Trixi.density(u, eqns::GramianMomentEquations1D3V) = u[1]


##########################################################
######################## Closure #########################
##########################################################
"""
    closure_moments(u, equations::GramianMomentEquations1D3V)

    All closing moments, i.e. the moments of total order M+1 appearing in the flux, in the
    order given by `equations._closure_index`.

    For each direction n(θ_i, φ_i) the state is rotated onto the directional moments
    ∫ f (v·n_i)^k dv, k = 0…M, the 1D Gramian closure supplies ∫ f (v·n_i)^(M+1) dv, and
    the resulting `N_shell` values are inverted onto the order-(M+1) moments. Only the
    α ≥ 1 prefix that the flux needs is returned.
"""
function closure_moments(u, equations::GramianMomentEquations1D3V{Mp1, N, MC, NC, NS}) where {Mp1, N, MC, NC, NS}
    T = eltype(u)
    directional = MVector{NS, T}(undef)
    @inbounds for i in 1:NS
        # Directional moments ∫ f (v·n_i)^k dv for k = 0…M …
        # ? Maybe the quicker way to use the definition from above:  u_rot = matvec(equations._P[i], u, Val(MC))
        u_rot = SMatrix{MC, length(u)}(equations._P[i]) * u
        # … and the 1D closure value ∫ f (v·n_i)^(M+1) dv.
        directional[i] = closure(u_rot, equations)
    end
    # ? return matvec(equations._T_closure, SVector(directional), Val(NC))
    return SMatrix{NC, NS}(equations._T_closure) * SVector(directional)
end

"""
    Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}

    Flux function: F(U_(αβγ)) = U_(α+1,βγ), where the moments of total order M+1 are
    obtained from the closure relation.
"""
function Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    closing_moments = closure_moments(u, equations)
    flux_index = equations._flux_index
    return SVector(ntuple(Val(Mp1)) do i
        j = flux_index[i]
        j > 0 ? u[j] : closing_moments[-j]
    end)
end


##########################################################
##################### Source terms #######################
##########################################################

function zero_source(u, x, t, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    return SVector(ntuple(i -> zero(eltype(u)), Val(Mp1)))
end

"""
    relaxation_source(u, x, t, equations::GramianMomentEquations1D3V)

    BGK relaxation towards the local Maxwellian, RHS = -1/Kn (u - u_eq).

    The equilibrium is built in central-moment space, where it is diagonal:
    P_eq_(αβγ) = ρ (α-1)!! (β-1)!! (γ-1)!! θ^((α+β+γ)/2) if α, β, γ are all even, else 0.
"""
function relaxation_source(u, x, t, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    T = eltype(u)
    prim = cons2prim(u, equations)
    ρ = prim[1]
    θ = temperature(prim, equations)

    p_eq = MVector{Mp1, T}(undef)
    p_eq[1] = ρ
    p_eq[2] = prim[2]; p_eq[3] = prim[3]; p_eq[4] = prim[4]   # velocity components
    @inbounds for idx in 1:Mp1
        i, j, k = equations._pow[idx]
        if i + j + k <= 1 # handle slab-geometry
            p_eq[idx] = prim[idx]   # ρ and the velocity components are shared with u_eq
        elseif isodd(i) || isodd(j) || isodd(k)
            p_eq[idx] = zero(T)
        else
            # (i-1)!! (j-1)!! (k-1)!! θ^((i+j+k)/2); the exponent is exact since i+j+k is even
            p_eq[idx] = ρ * double_factorial(i - 1) * double_factorial(j - 1) *
                        double_factorial(k - 1) * θ^((i + j + k) ÷ 2)
        end
    end

    u_eq = prim2cons(SVector(p_eq), equations)
    return SVector(ntuple(i -> -equations.inv_Kn * (u[i] - u_eq[i]), Val(Mp1)))
end


##########################################################
############### Conservative <-> primitive ###############
##########################################################

"""
    moment_position(equations, i, j, k)

    Position of U_(ijk) in the state vector, or 0 if that moment is not evolved.
"""
@inline function moment_position(equations::GramianMomentEquations1D3V, i::Int, j::Int, k::Int)
    (i < 0 || j < 0 || k < 0) && return 0
    (i > equations.M || j > equations.M || k > equations.M) && return 0
    return @inbounds equations._pos[i + 1, j + 1, k + 1]
end

"""
    moment_val(u, i, j, k, equations)

    U_(ijk) from the state vector `u`, or zero if that moment is not tracked.
"""
@inline function moment_val(u, i::Int, j::Int, k::Int, equations::GramianMomentEquations1D3V)
    p = moment_position(equations, i, j, k)
    return p == 0 ? zero(eltype(u)) : @inbounds u[p]
end

"""
    central_moment_val(prim, i, j, k, equations)

    P_(ijk) from a primitive vector. Slots 2:4 of `prim` hold the velocity components
    rather than the (identically vanishing) first-order central moments, so those are
    reported as zero, and P_(000) = ρ.
"""
@inline function central_moment_val(prim, i::Int, j::Int, k::Int, equations::GramianMomentEquations1D3V)
    s = i + j + k
    s == 0 && return @inbounds prim[1]           # P_(000) = ρ
    s == 1 && return zero(eltype(prim))          # P_(100) = P_(010) = P_(001) = 0
    p = moment_position(equations, i, j, k)
    return p == 0 ? zero(eltype(prim)) : @inbounds prim[p]
end

"""
    moment_cons2prim(u_cons, equations::GramianMomentEquations1D3V)

    Raw moments -> (ρ, v_x, v_y, v_z, central moments P_(αβγ)).

    P_(αβγ) = Σ_{l≤α, m≤β, p≤γ} C(α,l) C(β,m) C(γ,p) (-v_x)^(α-l) (-v_y)^(β-m) (-v_z)^(γ-p) U_(lmp)

    Unlike the slab system, the shift is performed in all three velocity directions.
"""
function moment_cons2prim(u_cons, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    T = eltype(u_cons)
    prim = MVector{Mp1, T}(undef)

    ρ = u_cons[1]
    # A velocity moment that is not evolved (slab mode) vanishes by symmetry.
    v_x = moment_val(u_cons, 1, 0, 0, equations) / ρ
    v_y = moment_val(u_cons, 0, 1, 0, equations) / ρ
    v_z = moment_val(u_cons, 0, 0, 1, equations) / ρ

    @inbounds for idx in 1:Mp1
        i, j, k = equations._pow[idx]
        s = i + j + k
        if s == 0
            prim[idx] = ρ
        elseif s == 1
            # the slots of the order-1 block hold the velocity components
            prim[idx] = (i == 1 ? v_x : (j == 1 ? v_y : v_z))
        else
            val = zero(T)
            for l in 0:i, m in 0:j, p in 0:k
                coef = binomial(i, l) * binomial(j, m) * binomial(k, p)
                val += coef * (-v_x)^(i - l) * (-v_y)^(j - m) * (-v_z)^(k - p) *
                    moment_val(u_cons, l, m, p, equations)
            end
            prim[idx] = val
        end
    end

    return SVector(prim)
end

"""
    moment_prim2cons(u_prim, equations::GramianMomentEquations1D3V)

    Inverse of [`moment_cons2prim`](@ref).
"""
function moment_prim2cons(u_prim, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    T = eltype(u_prim)
    cons = MVector{Mp1, T}(undef)

    ρ = u_prim[1]
    # In prim space the order-1 slots hold the velocities themselves; a component whose
    # moment is not evolved (slab mode) is zero by symmetry.
    p_y = moment_position(equations, 0, 1, 0)
    p_z = moment_position(equations, 0, 0, 1)
    v_x = u_prim[2]
    v_y = p_y == 0 ? zero(T) : @inbounds u_prim[p_y]
    v_z = p_z == 0 ? zero(T) : @inbounds u_prim[p_z]

    @inbounds for idx in 1:Mp1
        i, j, k = equations._pow[idx]
        s = i + j + k
        if s == 0
            cons[idx] = ρ
        elseif s == 1
            cons[idx] = ρ * (i == 1 ? v_x : (j == 1 ? v_y : v_z))
        else
            val = zero(T)
            for l in 0:i, m in 0:j, p in 0:k
                coef = binomial(i, l) * binomial(j, m) * binomial(k, p)
                val += coef * v_x^(i - l) * v_y^(j - m) * v_z^(k - p) *
                    central_moment_val(u_prim, l, m, p, equations)
            end
            cons[idx] = val
        end
    end

    return SVector(cons)
end

# Link to Trixi
Trixi.cons2prim(u, eqns::GramianMomentEquations1D3V) = moment_cons2prim(u, eqns)
Trixi.prim2cons(u, eqns::GramianMomentEquations1D3V) = moment_prim2cons(u, eqns)

# Convert conservative variables to entropy (necessary dummy)
Trixi.cons2entropy(u, equations::GramianMomentEquations1D3V) = u

##########################################################
#################### Maximum speeds ######################
##########################################################
const SPEED_FACTOR_1D3V = 5.0

"""
    temperature(prim, equations::GramianMomentEquations1D3V)

    θ = (P_(200) + P_(020) + P_(002)) / (3ρ) from a *primitive* (central-moment) vector.
"""
# ToDo: Check if this are not the conservative ones?
@inline function temperature(prim, equations::GramianMomentEquations1D3V)
    i200, i020, i002 = equations._pressure_index
    return @inbounds (prim[i200] + prim[i020] + prim[i002]) / (3 * prim[1])
end

"""
    λ_max = |v_x| + γ √θ, with γ = 5
"""
@inline function max_abs_speed_1d3v(u, equations::GramianMomentEquations1D3V)
    prim = cons2prim(u, equations)
    θ = temperature(prim, equations)
    return abs(prim[2]) + SPEED_FACTOR_1D3V * sqrt(max(θ, zero(θ)))
end

function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3V)
    return max(max_abs_speed_1d3v(u_l, equations), max_abs_speed_1d3v(u_r, equations))
end

function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D3V)
    return (max_abs_speed_1d3v(u, equations),)
end


##########################################################
################## Initial conditions ####################
##########################################################

"""
    Shock tube initial conditions with left and right distribution functions in 1D-3V.

    The moments are computed for the full multi-index list and then gathered into the
    ordering of `equations._U_t_index`.
"""
struct InitialConditionsShockTube1D3V{N}
    left::SVector{N, Float64}
    right::SVector{N, Float64}
end

function InitialConditionsShockTube1D3V(f_left, f_right, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    M = equations.M
    if equations.slab_geometry
        @assert f_left.v[2] == f_left.v[3] == 0 && f_right.v[2] == f_right.v[3] == 0 "slab geometry requires v_y = v_z = 0 in the initial data."
    end
    full = multi_index_list(M)
    position_in_full = Dict(ix => i for (i, ix) in enumerate(full))
    gather = [position_in_full[ix] for ix in equations._U_t_index]

    left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)[gather]
    right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)[gather]
    @assert length(left) == Mp1 == equations.N_equations

    return InitialConditionsShockTube1D3V{Mp1}(
        SVector{Mp1, Float64}(left),
        SVector{Mp1, Float64}(right)
    )
end

# Convenience method keeping the explicit `M` argument; it must agree with the equations.
function InitialConditionsShockTube1D3V(f_left, f_right, M::Integer, equations::GramianMomentEquations1D3V)
    @assert M == equations.M "initial condition built for M = $M but the equations use M = $(equations.M)."
    return InitialConditionsShockTube1D3V(f_left, f_right, equations)
end

function (ic::InitialConditionsShockTube1D3V)(coords, t, equations::GramianMomentEquations1D3V)
    return coords[1] < 0.0 ? ic.left : ic.right
end
