
"""
    Implementation of the Vlasov-Poisson system using Gramian moment equations in 1D1V.

    The electric field is computed from the density once per time step in a `DiscreteCallback`
    (`vlasov_poisson_callback`) by solving ∂E/∂x = ρ - ⟨ρ⟩ with periodic boundary conditions in
    Fourier space. It enters the moment equations through the source term `VlasovPoissonSource`.
    Buffers, FFT plans, the field interpolant and the energy history live in a `VlasovPoissonCache`
    that is created once per run and shared by the callback and the source term, so nothing is
    allocated inside the time loop.
"""

"""
    VlasovPoissonCache(mesh::TreeMesh{1}, solver::DGSEM)

    Preallocated workspace shared by `VlasovPoissonSource` and `vlasov_poisson_callback`.
    Requires a `LobattoLegendreBasis` of `polydeg` 1 or 2, see `check_uniform_vertex_grid`.

    # Fields:
    - `x::Vector{Float64}`: Coordinates of all `n+1` mesh vertices, `x[end] = x[1] + Lx`
    - `ρ::Vector{Float64}`: Density at the `n` periodic vertices (vertex 1 holds the periodic average)
    - `ρk::Vector{ComplexF64}`: Real-FFT spectrum of `ρ`, `n÷2+1` modes
    - `kx::Vector{Float64}`: Wavenumbers `2π/Lx * (0:n÷2)` matching `ρk`
    - `E_period::Vector{Float64}`: Inverse-FFT output, `n` samples
    - `E::Vector{Float64}`: Electric field at all `n+1` vertices, `E[end] = E[1]`
    - `E_itp`: Linear interpolant of `E` over `x`, re-wrapped once per step
    - `plan`, `iplan`: Precomputed forward and inverse real-FFT plans
    - `E_L2::Vector{Float64}`: L2 norm of the electric field at every step
    - `times::Vector{Float64}`: Times at which `E_L2` was recorded
"""
mutable struct VlasovPoissonCache{Itp, P, IP}
    x::Vector{Float64}
    ρ::Vector{Float64}
    ρk::Vector{ComplexF64}
    kx::Vector{Float64}
    E_period::Vector{Float64}
    E::Vector{Float64}
    E_itp::Itp
    plan::P
    iplan::IP
    E_L2::Vector{Float64}
    times::Vector{Float64}
end

"""
    check_uniform_vertex_grid(solver::DGSEM)

    The Poisson solve takes an FFT of the density sampled at the DG nodes, which only works for
    equispaced samples. Lobatto-Legendre nodes are equispaced for `polydeg` 1 (cell boundaries)
    and 2 (boundaries plus midpoint); from `polydeg` 3 on the interior nodes cluster towards the
    cell boundaries and the FFT would see a distorted grid (field error of 0.5 % at 64 cells,
    decaying only linearly with the mesh width, see docs/plans/vlasov_poisson_nonuniform_nodes.md). Throws an `ArgumentError` for any other basis or degree.
"""
function check_uniform_vertex_grid(solver::DGSEM)
    solver.basis isa LobattoLegendreBasis || throw(ArgumentError(
        "VlasovPoissonCache needs a LobattoLegendreBasis (nodes on the cell boundaries), got $(typeof(solver.basis))."))
    p = Trixi.polydeg(solver)
    p in (1, 2) || throw(ArgumentError(
        "VlasovPoissonCache needs polydeg 1 or 2: the FFT Poisson solve requires equispaced nodes, " *
        "and Lobatto-Legendre nodes are only equispaced up to polydeg 2 (got polydeg = $p)."))
    return nothing
end

# Global instance
function VlasovPoissonCache(mesh::TreeMesh{1}, solver::DGSEM)
    check_uniform_vertex_grid(solver)
    x = vertex_coordinates(mesh, solver)
    n = length(x) - 1                   # samples per period; x[end] is the periodic image of x[1]
    Lx = x[end] - x[1]
    ρ = zeros(n)
    ρk = zeros(ComplexF64, n ÷ 2 + 1)
    kx = collect((2π / Lx) .* (0:n÷2))
    E_period = zeros(n)
    E = zeros(n + 1)
    return VlasovPoissonCache(x, ρ, ρk, kx, E_period, E, field_interpolant(x, E), plan_rfft(ρ), plan_irfft(ρk, n), Float64[], Float64[])
end

"""
    field_interpolant(x, E)

    Linear interpolant of the electric field. `Interpolations` copies `E`, so this is rebuilt once
    per step in `solve_poisson!` (about 3 KB). Every DG node lies inside `x`; the extrapolation
    only guards against round-off at the domain ends.
"""
field_interpolant(x, E) = linear_interpolation(x, E, extrapolation_bc = Interpolations.Line())

"""
    vertex_coordinates(mesh::TreeMesh{1}, solver::DGSEM)

    Coordinates of all DG nodes of the 1D mesh from left to right, with nodes shared by two cells
    listed once (same enumeration as `collect1dTreeArrays`).
"""
function vertex_coordinates(mesh::TreeMesh{1}, solver::DGSEM)
    tree = mesh.tree
    nodes = solver.basis.nodes
    nn = length(nodes)
    cells = Trixi.leaf_cells(tree)
    x = zeros((nn - 1) * length(cells) + 1)
    for (i_cell, cell) in enumerate(cells)
        center = Trixi.cell_coordinates(tree, cell)[1]
        half_length = 0.5 * Trixi.length_at_cell(tree, cell)
        for i in 1:nn
            x[(i_cell - 1) * (nn - 1) + i] = center + half_length * nodes[i]
        end
    end
    return x
end

"""
    gather_density!(ρ, u_ode, semi)

    Density at the mesh vertices from the DG solution. Nodes shared by two cells are averaged,
    including across the periodic boundary, so `ρ` holds exactly one period.
"""
function gather_density!(ρ, u_ode, semi)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
    nn = Trixi.nnodes(solver)
    fill!(ρ, 0.0)
    for element in Trixi.eachelement(solver, cache), i in 1:nn
        vertex = mod1((element - 1) * (nn - 1) + i, length(ρ))  # last vertex wraps onto vertex 1
        weight = (i == 1 || i == nn) ? 0.5 : 1.0                 # cell boundary nodes are shared
        ρ[vertex] += weight * u[1, i, element]
    end
    return ρ
end

"""
    solve_poisson!(c::VlasovPoissonCache)

    Solve ∂E/∂x = ρ - ⟨ρ⟩ with periodic boundary conditions in Fourier space, in place.
    Afterwards `c.E` holds the field at all vertices and `c.E_itp` interpolates it.
"""
function solve_poisson!(c::VlasovPoissonCache)
    c.ρ .-= mean(c.ρ)                       # neutralizing background
    mul!(c.ρk, c.plan, c.ρ)                 # ρ̂(k)
    c.ρk[1] = 0                             # k = 0 mode: E has zero mean
    @inbounds for j in 2:length(c.ρk)
        c.ρk[j] /= im * c.kx[j]             # Ê(k) = ρ̂(k) / (i k)
    end
    mul!(c.E_period, c.iplan, c.ρk)         # back to physical space
    c.E[1:end-1] .= c.E_period
    c.E[end] = c.E[1]                       # periodic image of the first vertex
    c.E_itp = field_interpolant(c.x, c.E)
    return c
end

"""
    VlasovPoissonSource(cache::VlasovPoissonCache)

    Source term of the moment equations due to the electric field, `source_k = (k-1) E u_{k-1}`
    for `k ≥ 2` and no source for mass conservation. Pass it to `SemidiscretizationHyperbolic`
    as `source_terms`.
"""
struct VlasovPoissonSource{C <: VlasovPoissonCache}
    cache::C
end

function (s::VlasovPoissonSource)(u, x, t, equations::GramianMomentEquations1D{Mp1}) where {Mp1}
    E_local = s.cache.E_itp(x[1])
    # positive sign: the term sits on the left-hand side with a minus
    return SVector{Mp1}(ntuple(k -> k == 1 ? 0.0 : (k - 1) * E_local * u[k - 1], Val(Mp1)))
end

"""
    update_electric_field!(c::VlasovPoissonCache, integrator)

    Recompute the electric field from the current density and record `‖E‖_L2`.
"""
function update_electric_field!(c::VlasovPoissonCache, integrator)
    gather_density!(c.ρ, integrator.u, integrator.p)
    solve_poisson!(c)
    push!(c.E_L2, sqrt(trapz(c.x, c.E .^ 2)))   # ‖E(t,⋅)‖_L2 via the trapezoidal rule
    push!(c.times, integrator.t)
    return nothing
end

"""
    vlasov_poisson_callback(cache::VlasovPoissonCache)

    `DiscreteCallback` that recomputes the electric field after every step and once before the
    first one. The field history is available afterwards in `cache.E_L2` and `cache.times`.
"""
function vlasov_poisson_callback(cache::VlasovPoissonCache)
    update_field!(integrator) = update_electric_field!(cache, integrator)
    return DiscreteCallback(
        (u, t, integrator) -> true,                 # fire after every step
        update_field!;
        save_positions = (false, false),            # u is not modified, nothing extra to save
        initialize = (c, u, t, integrator) -> update_field!(integrator),
    )
end