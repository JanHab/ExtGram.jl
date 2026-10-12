import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ExtGram, LinearAlgebra

# moments_counterexample.jl
# Raw moments u_0..u_{2n-1} of the odd-M extended Gramian counterexample
# (frame a_{n-1} = 0, b_{n-1} = 1), exact in Rational{BigInt}.

"""
    moments_from_recurrence(a, b; u0 = 1)

a = (a_0, …, a_{n-1}), b = (b_1, …, b_{n-1}) — Jacobi recurrence coefficients
p_{k+1}(c) = (c - a_k) p_k(c) - b_k p_{k-1}(c).  Returns u_0 … u_{2n-1} with
u_k = u0 * e_1' J^k e_1, where J is the (non-symmetric, similar) Jacobi matrix
with J[k+1,k] = b_k and J[k,k+1] = 1 — same spectrum, no square roots.
"""
function moments_from_recurrence(a::AbstractVector{T}, b::AbstractVector{T}; u0::T = one(T)) where T
    n = length(a)
    @assert length(b) == n - 1
    J = zeros(T, n, n)
    for k in 1:n
        J[k, k] = a[k]
        if k ≥ 2
            J[k, k-1] = b[k-1]      # b_{k-1} on the sub-diagonal
            J[k-1, k] = one(T)      # 1 on the super-diagonal
        end
    end
    v = zeros(T, n); v[1] = one(T)
    u = Vector{T}(undef, 2n)
    for k in 1:2n
        u[k] = u0 * v[1]            # u_{k-1}
        v = J * v
    end
    return u
end

# ---------------------------------------------------------------- M = 21 counterexample
M   = 21 
Kn  = 1.0 # Doesn't really matter, as only the moments are computed here, not the PDE solved
closure = "ExtGram"
n   = (M + 1) ÷ 2
σ   = 1//4
ε   = 1//1000
a   = Rational{BigInt}[fill(-σ, n-1); 0]   # a_0 … a_9 = -1/4,  a_10 = 0
b   = Rational{BigInt}[fill(ε, n-2); 1]    # b_1 … b_9 = 1/1000, b_10 = 1

u = moments_from_recurrence(a, b)          # u_0 … u_21, u_0 = 1

for (k, uk) in enumerate(u)
    println("u_", k-1, " = ", uk, "   ≈ ", Float64(uk))
end

equations = GramianMomentEquations1D(M, Kn, closure)
u_closure = ExtGram.closure(u, equations)  # compute closure u_{22}
u_extended = vcat(u, u_closure)              # u_0 … u_22
G_n = ExtGram.gramian(u_extended, n)
# Check that moments are realizable
@assert all(eigen(G_n).values .> 0) "Gramian is not positive definite"

e_vals_flux = eigen(ExtGram.flux_jacobian(u, equations)).values
maximum_imaginary = maximum(abs.(imag.(e_vals_flux)))
if maximum_imaginary > 0
    @warn "Flux Jacobian has complex eigenvalues, hyperbolicity is lost."
    println("Maximum imaginary part of eigenvalues of flux Jacobian: ", maximum_imaginary)
else
    println("All eigenvalues of flux Jacobian are real.")
end





M   = 19 
Kn  = 1.0 # Doesn't really matter, as only the moments are computed here, not the PDE solved
closure = "ExtGram"
n   = (M + 1) ÷ 2
σ   = 1//4
ε   = 1//1000
a   = Rational{BigInt}[fill(-σ, n-1); 0]   # a_0 … a_9 = -1/4,  a_10 = 0
b   = Rational{BigInt}[fill(ε, n-2); 1]    # b_1 … b_9 = 1/1000, b_10 = 1

u = moments_from_recurrence(a, b)          # u_0 … u_21, u_0 = 1

for (k, uk) in enumerate(u)
    println("u_", k-1, " = ", uk, "   ≈ ", Float64(uk))
end

equations = GramianMomentEquations1D(M, Kn, closure)
u_closure = ExtGram.closure(u, equations)  # compute closure u_{22}
u_extended = vcat(u, u_closure)              # u_0 … u_22
G_n = ExtGram.gramian(u_extended, n)
# Check that moments are realizable
@assert all(eigen(G_n).values .> 0) "Gramian is not positive definite"

e_vals_flux = eigen(ExtGram.flux_jacobian(u, equations)).values
maximum_imaginary = maximum(abs.(imag.(e_vals_flux)))
if maximum_imaginary > 0
    @warn "Flux Jacobian has complex eigenvalues, hyperbolicity is lost."
    println("Maximum imaginary part of eigenvalues of flux Jacobian: ", maximum_imaginary)
else
    println("All eigenvalues of flux Jacobian are real.")
end