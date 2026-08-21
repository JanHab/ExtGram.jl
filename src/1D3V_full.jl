
function U_t_index(nmax::Integer, slab_geometry::Bool=false)
    if slab_geometry
        # ToDo: Rename get_valid_indices to something more descriptive, e.g., get_slab_indices
        return mainmomindex(nmax)[ExtGram.get_valid_indices(nmax)]
    else
        return mainmomindex(nmax)
    end
end

function U_x_index(nmax::Integer, slab_geometry::Bool=false)
    U_t = U_t_index(nmax, slab_geometry)
    return [U_t[i] .+ [1, 0, 0] for i in 1:length(U_t)]
end



struct GramianMomentEquations1D3V{Mp1, N, NC, RealT <: Real} <: GramianMomentEquations{Mp1, N, RealT}
    inv_Kn::RealT # 1/Kn (Kn = Knudsen number) | used by the relaxation source term
    χ::RealT
    M::Int          # highest moment order
    n::Int
    N_equations::Int
    N_closure::Int  # number of moments the closure has to supply
    closure::Symbol
    _U_t_index::Vector{Vector{Int}}
    _U_x_index::Vector{Vector{Int}}
    # Flux lookup table, one entry per equation i:
    #   j > 0  ->  U_(_U_x_index[i]) is itself an evolved moment, found at u[j]
    #   j < 0  ->  U_(_U_x_index[i]) is a closing moment, found at closure_moments[-j]
    _flux_index::Vector{Int}
    # Multi-indices of the closing moments, in the order the closure must return them
    _closure_index::Vector{Vector{Int}}
    # Rotational matrix for the closure moments
    T_tilde_inv::MMatrix{NC, NC, RealT}
    # Angles for the rotational matrix
    theta::Vector{RealT}
    phi::Vector{RealT}

    function GramianMomentEquations1D3V(
        M::Integer, Knudsen::Real, closure::String="ExtGram"; 
        slab_geometry::Bool=false,
        theta::Vector{<:Real}=Float64[], phi::Vector{<:Real}=Float64[]
    )
        @assert M > 1
        @assert length(theta) == length(phi)
        if iseven(M)
            n = Int(M/2)
            χ = (n+1)/n
        else
            n = Int((M+1)/2)
            χ = (n+1)/(2n)
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
        _flux_index = Vector{Int}(undef, N_equations)
        _closure_index = Vector{Vector{Int}}()
        for i in 1:N_equations
            j = get(position_in_u, _U_x_index[i], 0)
            if j != 0
                _flux_index[i] = j
            else
                push!(_closure_index, _U_x_index[i])
                _flux_index[i] = -length(_closure_index)
            end
        end
        N_closure = length(_closure_index)

        # The rotational matrix
        RealT = typeof(Knudsen)
        T_tilde = zeros(MMatrix{N_closure, N_closure, RealT})
        for i in 1:length(theta)
            T_tilde[i, :] = ExtGram.tensor_transformation(M, theta[i], phi[i])[1, :]
        end

        # NOTE: the first type parameter is Trixi's NVARS (see the abstract type in
        # src/gramian_moment_equations.jl), so it must be N_equations, not M+1.
        new{N_equations, n, N_closure, typeof(Knudsen)}(
            inv(Knudsen), χ, M, n, N_equations, N_closure,
            closure_value, _U_t_index, _U_x_index,
            _flux_index, _closure_index, inv(T_tilde),
            convert(Vector{RealT}, theta),
            convert(Vector{RealT}, phi)
        )
    end
end

# # * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), equations::GramianMomentEquations1D3V) = ntuple(i->"U_($(equations._U_t_index[i]))", equations.N_equations)
Trixi.varnames(::typeof(cons2prim), equations::GramianMomentEquations1D3V) = ntuple(i->"U_($(equations._U_x_index[i]))", equations.N_equations)

# * basic variable definitions to be used as e.g. sc indicator variable
Trixi.density(u, eqns::GramianMomentEquations1D3V) = u[1]

"""
    closure_moments(u, equations::GramianMomentEquations1D3V{Mp1, N, NC}) where {Mp1, N, NC}

    Returns all closing moments, i.e. the moments of total order M+1 appearing in the flux,
    in the order given by `equations._closure_index`.

    ToDo: PLACEHOLDER — returns zeros so that the flux assembly can be tested.
"""
function closure_moments(u, equations::GramianMomentEquations1D3V{Mp1, N, NC}) where {Mp1, N, NC}
    U_tilde_full = zeros(MVector{length(equations._closure_index)})
    for i in 1:length(equations.theta)
        theta_ = equations.theta[i]
        phi_ = equations.phi[i]
        # println("theta: $theta_, phi: $phi_")

        u_rot = MVector{length(u)}(u)
        lower_index = 1
        for M_ = 0:equations.M
            # ToDo: Once compute t_transformations for all M_ in advance, we can just slice the precomputed matrices here.
            t_transformations = ExtGram.tensor_transformation(M_, theta_, phi_)
            upper_index = lower_index + size(t_transformations)[1] # it is a square matrix, so size(t_transformations)[1] == size(t_transformations)[2]
            u_rot[lower_index:upper_index-1] = t_transformations * u_rot[lower_index:upper_index-1]

            lower_index = upper_index
        end

        # U_{α, 0, 0} indizes
        U_alpha_0_0_indices = ExtGram.nidx(M)
        U_tilde_Mp1 = closure(u_rot[U_alpha_0_0_indices], equations)
        U_tilde_full[i] = U_tilde_Mp1
    end

    U_Mp1 = equations.T_tilde_inv * U_tilde_full

    return U_Mp1
end

"""
    Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}

    Flux function: F(U_{α, β, γ}) = U_{α+1, β, γ}
    where α, β, γ are the indices of the moments in 3D.
    The last moments are obtained from the closure relation.
"""
function Trixi.flux(u, orientation::Integer, equations::GramianMomentEquations1D3V{Mp1}) where {Mp1}
    closing_moments = closure_moments(u, equations)
    flux_index = equations._flux_index
    return SVector(ntuple(Val(Mp1)) do i
        j = flux_index[i]
        j > 0 ? u[j] : closing_moments[-j]
    end)
end

function zero_source(u, x, t, equations::GramianMomentEquations1D3V)
    return SVector{equations.N_equations}(ntuple(i->0.0, equations.N_equations))
end

"""
    get_moment_val(u, i, j, k, equations::GramianMomentEquations1D3V)

Helper to retrieve U_{ijk} from the state vector `u` based on the provided 
indices_list. Returns 0.0 if the moment is not in the system.
"""
@inline function get_moment_val(u, i_req, j_req, k_req, equations::GramianMomentEquations1D3V)
    # Search for the tuple in the index list (constant folding should optimize this)
    for (idx, (i, j, k)) in enumerate(equations._U_t_index)
        if i == i_req && j == j_req && k == k_req
            return u[idx]
        end
    end
    return 0.0 # Return 0 if moment not tracked (or handle error)
end

function moment_cons2prim(u_cons, equations::GramianMomentEquations1D3V)
    # Type and length of array
    T = eltype(u_cons)
    N = length(equations._U_t_index)

    # Initialize Mutable Vector
    prim = MVector{N, T}(undef)

    # Explicitly set the first two moments (Avoiding the if-check)
    prim[1] = u_cons[1]; prim[2] = u_cons[2] / prim[1] # rho, v

    # Loop for the remaining moments (Indices 3 to N)
    v = prim[2]
    for idx in 3:N
        i, j, k = equations._U_t_index[idx]
        val = zero(T)
        
        # Calculate P_{ijk} = sum Binomial * (-v)^... * U
        for m in 0:i
            u_mjk = get_moment_val(u_cons, m, j, k, equations)
            val += binomial(i, m) * (-v)^(i-m) * u_mjk
        end
        
        prim[idx] = val
    end

    return SVector(prim)
end

# Link to Trixi
Trixi.cons2prim(u, eqns::GramianMomentEquations1D3V) = moment_cons2prim(u, eqns)

# This remains the same, but now calls the faster constant version above.
"""
    λ_max = v_x ± γ √(θ), with γ=5
"""
# ToDo: Rewrite if reduce to 1D-2V or 1D-3V slab-geometry
function Trixi.max_abs_speed_naive(u_l, u_r, orientation::Integer, equations::GramianMomentEquations1D3V)
    ρ_l = u_l[1]; ρ_r = u_r[1]

    # Calculate Central Moments P (needed for Temperature)
    p_prim_l = cons2prim(u_l, equations); p_prim_r = cons2prim(u_r, equations);
    
    P_200_l = p_prim_l[5]; P_200_r = p_prim_r[5] # P_xx
    P_020_l = p_prim_l[8]; P_020_r = p_prim_r[8] # P_yy
    P_002_l = p_prim_l[10]; P_002_r = p_prim_r[10] # P_zz

    v_kl = p_prim_l[2]; v_kr = p_prim_r[2] # velocity in x-direction
    
    # Temperature theta = (U_200 + U_020 + U_002) / (3 * rho)
    θ_l = 1 / (3 * ρ_l) * (P_200_l + P_020_l + P_002_l - ρ_l * v_kl^2)
    θ_r = 1 / (3 * ρ_r) * (P_200_r + P_020_r + P_002_r - ρ_r * v_kr^2)
    
    γ = 5.0
    λ_l = ρ_l + γ * sqrt(θ_l)
    λ_r = ρ_r + γ * sqrt(θ_r)
    return max(λ_l, λ_r)
end

# ? How can we use the one from above?
function Trixi.max_abs_speeds(u, equations::GramianMomentEquations1D3V)
    # For the moment hard coded, not ideal
    # return 10.0 

    ρ = u[1]

    # Calculate Central Moments P (needed for Temperature)
    p_prim = cons2prim(u, equations)
    
    P_200 = p_prim[5] # P_xx
    P_020 = p_prim[8] # P_yy
    P_002 = p_prim[10] # P_zz

    v_k = p_prim[2] # velocity in x-direction
    
    # Temperature theta = (U_200 + U_020 + U_002) / (3 * rho)
    θ = 1 / (3 * ρ) * (P_200 + P_020 + P_002 - ρ * v_k^2)
    
    γ = 5.0
    λ = ρ_ + γ * sqrt(θ)

    return λ
end





# * Some placeholder for innitialization

function raw_moment_1d(k::Integer, u::Real, σ::Real)
    m_prev, m = one(u), u          # m_0, m_1
    k == 0 && return m_prev
    for j in 1:(k - 1)
        m_prev, m = m, u * m + j * σ^2 * m_prev
    end
    return m
end

# maxwellian_moment(idx, ρ, u, σ) =
#     ρ * prod(raw_moment_1d(idx[d], u[d], σ[d]) for d in 1:3)

# const ρ_test = 1.7
# const u_test = (0.35, -0.8, 0.15)
# const σ_test = (1.1, 0.7, 1.3)

# maxwellian_moments(indices) = [maxwellian_moment(i, ρ_test, u_test, σ_test) for i in indices]

# shift_x(idx) = idx .+ [1, 0, 0]

struct InitialConditionsShockTube1D3V{N}
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsShockTube1D3V(f_left, f_right, M, eqns::GramianMomentEquations1D3V{Mp1}) where {Mp1}
        # Compute convective moments
        left = convective_moments_1D3D(M, f_left.ρ, f_left.v, f_left.θ)
        right = convective_moments_1D3D(M, f_right.ρ, f_right.v, f_right.θ)

        return new{length(left)}(left, right)
    end
end

function (ic::InitialConditionsShockTube1D3V)(coords, t, equations::GramianMomentEquations1D3V)
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end