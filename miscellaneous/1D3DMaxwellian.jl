using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using Printf
using Plots

# 1D-space / 3D-velocity Maxwellian utilities
# Functions here compute the Maxwellian value and its velocity moments
# using tensor-product Gauss–Hermite quadrature.

const TWO_PI = 2π

"""
    maxwellian_3d(vx, vy, vz; rho=1.0, u=(ux,uy,uz), theta=1.0)

Return the Maxwellian PDF value at velocity (vx,vy,vz):
  f(v) = rho / (2π θ)^(3/2) * exp(-|v - u|^2 / (2 θ))
Parameters are keyword args: `rho`, `u` (3-tuple), and `theta` (temperature).
"""
function maxwellian_3d(vx::Real, vy::Real, vz::Real; rho::Real=1.0, u::NTuple{3,Real}=(0.0,0.0,0.0), theta::Real=1.0)
    ux, uy, uz = u
    dx2 = (vx - ux)^2
    dy2 = (vy - uy)^2
    dz2 = (vz - uz)^2
    pref = rho / (TWO_PI * theta)^(3/2)
    return pref * exp(-(dx2 + dy2 + dz2) / (2 * theta))
end

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
    multi_index_list(max_degree)

Return a flattened list of multi-indices for all shells 0..max_degree.
"""
function multi_index_list(max_degree::Integer)
    return vcat([index3D(n) for n in 0:max_degree]...)
end

"""
    velocity_moments_3d(max_degree; rho=1.0, u=(0,0,0), theta=1.0, q=8)

Compute the velocity moments up to `max_degree` (all multi-indices with i+j+k <= max_degree)
using tensor-product Gauss–Hermite with `q` nodes per dimension. Returns a Vector{Float64}
with the same ordering as `multi_index_list(max_degree)`.
"""
function velocity_moments_3d(max_degree::Integer; rho::Real=1.0, u::NTuple{3,Real}=(0.0,0.0,0.0), theta::Real=1.0, q::Integer=8)
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
    idxs = multi_index_list(max_degree)
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

# ---------------- Example usage ----------------
println("Building Maxwellian moments example")
rho = 1.0
u = (0.2, 0.0, 0.0)
theta = 1.0
maxdeg = 4 # maximum total degree of moments to compute
q = 10 # number of quadrature points per velocity dimension

# --- Visualization: 1D slice (vx) and 2D heatmap (vx,vy) at vz=0 ---
function plot_maxwellian_slice(rho, u::NTuple{3,Real}, theta; nv=401, range_factor=4, vy=0.0, vz=0.0, savepath="miscellaneous/maxwellian_slice.png")
    ux = u[1]
    r = range(ux - range_factor*sqrt(theta), ux + range_factor*sqrt(theta), length=nv)
    vals = [maxwellian_3d(x, vy, vz; rho=rho, u=u, theta=theta) for x in r]
    plt = plot(r, vals, xlabel="vx", ylabel="f(vx,vy=0,vz=0)", title="Maxwellian slice (vx)", legend=false)
    # savefig(plt, savepath)
    display(plt)
    println("Saved 1D Maxwellian slice to ", savepath)
end

function plot_maxwellian_heatmap(rho, u::NTuple{3,Real}, theta; nv=201, range_factor=3, vz=0.0, savepath="miscellaneous/maxwellian_heatmap.png")
    ux, uy, _ = u
    rx = range(ux - range_factor*sqrt(theta), ux + range_factor*sqrt(theta), length=nv)
    ry = range(uy - range_factor*sqrt(theta), uy + range_factor*sqrt(theta), length=nv)
    F = zeros(Float64, nv, nv)
    for (i,x) in enumerate(rx), (j,y) in enumerate(ry)
        F[j,i] = maxwellian_3d(x, y, vz; rho=rho, u=u, theta=theta) # rows->y for plotting orientation
    end
    plt = heatmap(rx, ry, F, xlabel="vx", ylabel="vy", title="Maxwellian (vz=0)")
    # savefig(plt, savepath)
    display(plt)
    println("Saved 2D Maxwellian heatmap to ", savepath)
end

# Create visualizations before computing moments
plot_maxwellian_slice(rho, u, theta)
plot_maxwellian_heatmap(rho, u, theta)

# Compute moments after visualization
moms = velocity_moments_3d(maxdeg; rho=rho, u=u, theta=theta, q=q)
idxs = multi_index_list(maxdeg)
println("Computed moments up to degree ", maxdeg)
for (m,(i,j,k)) in zip(moms, idxs)
    @printf("(%d,%d,%d): %.6g\n", i, j, k, m)
end



# Example IC:
x = LinRange(-2.0, 2.0, 1_000)
ρ1 = 7.0
ρ2 = 1.0
u1 = (0.2, 0.0, 0.0)
u2 = (0.0, 0.0, 0.0)
theta1 = 1.0
theta2 = 1.0
maxdeg = 4 # maximum total degree of moments to compute
q = 10 # number of quadrature points per velocity dimension

moments = [
    x_i < 0.0 ? velocity_moments_3d(maxdeg; rho=ρ1, u=u1, theta=theta1, q=q) : velocity_moments_3d(maxdeg; rho=ρ2, u=u2, theta=theta2, q=q) for x_i in x]

# Plotting density moment (0,0,0)
plt = plot()
for i in 1:size(moments[1], 1)
    moment = [m[i] for m in moments]
    plot!(plt, x, moment, xlabel="x", ylabel="u", label="u$(i)")
end
display(plt)