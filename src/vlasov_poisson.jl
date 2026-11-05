
"""
    Implementation of the Vlasov-Poisson system using Gramian moment equations in 1D1V.

    The electric field is computed from the density of the state being differentiated, once per
    RHS evaluation (`rhs_vlasov_poisson!`, installed by `semidiscretize_vlasov_poisson`), by
    solving ∂E/∂x = ρ - ⟨ρ⟩ with periodic boundary conditions in Fourier space. It enters the
    moment equations through the source term `VlasovPoissonSource`, which reads the field at a
    DG node by direct lookup on the equispaced vertex grid (`field_at`). `vlasov_poisson_callback`
    only records the field energy ½∫E²dx of the accepted states. Buffers, FFT plans and the energy
    history live in a
    `VlasovPoissonCache` that is created once per run and shared by the RHS, the source term and
    the callback, so nothing is allocated inside the time loop.
"""

"""
    VlasovPoissonCache(mesh::TreeMesh{1}, solver::DGSEM)

    Preallocated workspace shared by `rhs_vlasov_poisson!`, `VlasovPoissonSource` and
    `vlasov_poisson_callback`. Requires a `LobattoLegendreBasis` of `polydeg` 1 or 2, see
    `check_uniform_vertex_grid`.

    # Fields:
    - `x::Vector{Float64}`: Coordinates of all `n+1` mesh vertices, `x[end] = x[1] + Lx`
    - `x₁::Float64`, `h::Float64`, `n::Int`: First vertex, vertex spacing and number of periodic
      vertices, used by the field lookup `field_at`
    - `ρ::Vector{Float64}`: Density at the `n` periodic vertices (vertex 1 holds the periodic average)
    - `ρk::Vector{ComplexF64}`: Real-FFT spectrum of `ρ`, `n÷2+1` modes
    - `kx::Vector{Float64}`: Wavenumbers `2π/Lx * (0:n÷2)` matching `ρk`
    - `E_period::Vector{Float64}`: Inverse-FFT output, `n` samples
    - `E::Vector{Float64}`: Electric field at all `n+1` vertices, `E[end] = E[1]`
    - `plan`, `iplan`: Precomputed forward and inverse real-FFT plans
    - `E_energy::Vector{Float64}`: Field energy ½∫E²dx at t = 0 and after every accepted step
    - `times::Vector{Float64}`: Times at which `E_energy` was recorded
"""
struct VlasovPoissonCache{P, IP}
    x::Vector{Float64}
    x₁::Float64
    h::Float64
    n::Int
    ρ::Vector{Float64}
    ρk::Vector{ComplexF64}
    kx::Vector{Float64}
    E_period::Vector{Float64}
    E::Vector{Float64}
    plan::P
    iplan::IP
    E_energy::Vector{Float64}
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
    return VlasovPoissonCache(x, x[1], Lx / n, n, ρ, ρk, kx, E_period, E, plan_rfft(ρ), plan_irfft(ρk, n), Float64[], Float64[])
end

"""
    field_at(c::VlasovPoissonCache, x)

    Electric field at the coordinate `x` by linear lookup on the equispaced vertex grid. Same value
    as a linear interpolant of `c.E` over `c.x` with linear extrapolation, without the copy of `E`
    that `Interpolations.jl` makes on every call. Every DG node coincides with a vertex, so the
    lookup is exact there; at `c.x[end]` it returns `E[1]`; outside the grid (round-off at the
    domain ends only) it extrapolates along the end segments. Allocation-free.
"""
@inline function field_at(c::VlasovPoissonCache, x)
    σ = (x - c.x₁) / c.h
    i = clamp(floor(Int, σ), 0, c.n - 1)                 # segment [x[i+1], x[i+2]]
    @inbounds return c.E[i + 1] + (σ - i) * (c.E[i + 2] - c.E[i + 1])
end

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
    Afterwards `c.E` holds the field at all vertices.
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
    E_local = field_at(s.cache, x[1])
    # positive sign: the term sits on the left-hand side with a minus
    return SVector{Mp1}(ntuple(k -> k == 1 ? 0.0 : (k - 1) * E_local * u[k - 1], Val(Mp1)))
end

"""
    update_electric_field!(c::VlasovPoissonCache, u_ode, semi)

    Gather the density of `u_ode` and solve for the electric field into `c`. No history, no
    allocation. Called once per RHS evaluation and once per accepted step by the callback.
"""
function update_electric_field!(c::VlasovPoissonCache, u_ode, semi)
    Trixi.@trixi_timeit Trixi.timer() "Poisson solve" begin   # shows up in SummaryCallback
        gather_density!(c.ρ, u_ode, semi)
        solve_poisson!(c)
    end
    return c
end

"""
    rhs_vlasov_poisson!(du_ode, u_ode, semi, t)

    ODE right-hand side: field of *this* state, then the Trixi RHS.
    Replaces `Trixi.rhs!` as the ODE function so that the field is recomputed for every RK
    stage, not once per accepted step. `vlasov_poisson_callback` only records the field energy.
"""
function rhs_vlasov_poisson!(du_ode, u_ode, semi, t)
    update_electric_field!(semi.source_terms.cache, u_ode, semi)   # cache lives in the source
    Trixi.rhs!(du_ode, u_ode, semi, t)
    return nothing
end

"""
    semidiscretize_vlasov_poisson(semi, tspan; kwargs...)

    Drop-in replacement for `semidiscretize` that installs the wrapped RHS.
"""
function semidiscretize_vlasov_poisson(semi, tspan; kwargs...)
    ode = semidiscretize(semi, tspan; kwargs...)
    # Mirror Trixi: in place and fully specialised. A plain `ODEProblem(f, u0, tspan, p)` is
    # AutoSpecialize and wraps `f` in a FunctionWrapper. `ode.p` rather than `semi`: with
    # `storage_type`/`real_type` Trixi replaces the semidiscretization.
    return ODEProblem{true, OrdinaryDiffEq.SciMLBase.FullSpecialize}(
        rhs_vlasov_poisson!, ode.u0, ode.tspan, ode.p)
end

"""
    record_field!(c::VlasovPoissonCache, integrator)

    Called by the callback after every accepted step: E of the accepted state, then the field
    energy ½∫E²dx (trapezoidal rule on the vertex grid, exact for the FFT field) and t.
"""
function record_field!(c::VlasovPoissonCache, integrator)
    update_electric_field!(c, integrator.u, integrator.p)
    push!(c.E_energy, 1/2 * trapz(c.x, c.E .^ 2))   # ½∫E²dx
    push!(c.times, integrator.t)
    u_modified!(integrator, false)   # u untouched: keep the FSAL stage, no extra RHS call
    return nothing
end

"""
    check_wrapped_rhs(integrator)

    Refuse the lagged coupling: the ODE must have been built by `semidiscretize_vlasov_poisson`.
"""
function check_wrapped_rhs(integrator)
    f = OrdinaryDiffEq.SciMLBase.unwrapped_f(integrator.f.f)   # sees through FunctionWrappers
    f === rhs_vlasov_poisson! || throw(ArgumentError(
        "vlasov_poisson_callback needs the ODE from semidiscretize_vlasov_poisson; " *
        "with plain semidiscretize every RK stage would use the field of the previous step."))
    return nothing
end

"""
    vlasov_poisson_callback(cache::VlasovPoissonCache)

    DiscreteCallback that recomputes the electric field and records its energy ½∫E²dx after every
    accepted step. The callback is installed in `solve` after `semidiscretize_vlasov_poisson`.
"""
function vlasov_poisson_callback(cache::VlasovPoissonCache)
    affect!(integrator) = record_field!(cache, integrator)
    return DiscreteCallback(
        (u, t, integrator) -> true,                 # fire after every accepted step
        affect!;
        save_positions = (false, false),
        initialize = (c, u, t, integrator) -> (check_wrapped_rhs(integrator); affect!(integrator)),
    )
end
