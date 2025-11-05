
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
