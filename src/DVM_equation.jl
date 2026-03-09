struct DVMEquations1D{N} <: Trixi.AbstractEquations{1, N}
    """
    Discrete Velocity Method (DVM) Equations in 1D velocity space with N discrete velocity points.

    ∂_t f_k + c_k ∂_x f_k = - (1/Kn) (f_k - f̃_k),  k = 1,...,N

    where f̃_k is the corrected Maxwellian distribution ensuring conservation of mass, momentum, and energy.
    Discretizes the velocity space using N discrete velocity points between c_l and c_u.

    # Fields
    - `N::Int`: Number of equations to be solved.
    - `c_l::Real`: Lower bound of the velocity domain.
    - `c_u::Real`: Upper bound of the velocity domain.
    - `c_vec::SVector{N, Float64}`: Velocity grid points.
    - `Kn::Real`: Knudsen number.
    """
    N::Int  # Number of equations to be solved
    c_l::Real   # Lower bound of the velocity domain
    c_u::Real   # Upper bound of the velocity domain
    c_vec::SVector{N, Float64}  # Velocity grid points
    Kn::Real    # Knudsen number
    function DVMEquations1D(N::Integer, c_l::Real, c_u::Real, Kn::Real)
        @assert N > 1 "Number of equations N must be positive"
        @assert c_l < c_u "Lower bound c_l must be smaller than upper bound c_u"
        c_vec = SVector{N, Float64}(Tuple(LinRange(c_l, c_u, N)))
        new{N}(N, c_l, c_u, c_vec, Kn)
    end
end
dc(equations::DVMEquations1D) = (equations.c_u - equations.c_l) / (equations.N - 1)

# * Conservative to Primitive variables and vice versa
Trixi.varnames(::typeof(cons2cons), ::DVMEquations1D{N}) where {N} = ntuple(i->"f^{$(i-1)}", N)

function Trixi.flux(f, orientation::Integer, equations::DVMEquations1D{N}) where {N}
    """
    Compute the flux for the DVM equations in 1D velocity space.

    f_i * c_i for i = 1,...,N
    """
    return SVector(
        ntuple(i -> f[i] * equations.c_vec[i], N)
    )
end

function relaxation_source(f, x, t, equations::DVMEquations1D{N}) where {N} 
    """
    Compute the DVM relaxation source term with conservation constraints.

    Relaxation is: (-1/Kn) * (f - f̃)

    where f̃ is the corrected Maxwellian distribution ensuring conservation of mass, momentum, and energy.
    ρ = ∫ f dc
    ρ v = ∫ c f dc
    ρ θ = ∫ (c - v)² f dc

    Steps:
    1. Compute macroscopic moments ρ, v, θ from f.
    2. Compute the Maxwellian distribution f_Maxwellian based on these moments.
    3. Set up the Lagrange-multiplier system to enforce conservation of mass, momentum, and energy.
    4. Solve for the correction Δf to f_Maxwellian.
    5. Return the relaxation source term.
    """
    ρ = trapz(equations.c_vec, f) # density
    v = trapz(equations.c_vec, f .* equations.c_vec) / ρ # velocity
    θ = trapz(equations.c_vec, f .* (equations.c_vec .- v).^ 2) / ρ # temperature

    ρ = max.(ρ, 1e-12)
    θ = max.(θ, 1e-12)
    f_Maxwellian = Maxwellian(ρ, v, θ).(equations.c_vec)
    return (-1.0 / equations.Kn) * (f .- f_Maxwellian)

end

# Zero source term (collisionless case)
zero_source(f, x, t, equations::DVMEquations1D{N}) where {N} = SVector{N}(ntuple(i->0.0, N))

function flux_jacobian(u, equations::DVMEquations1D{N}) where {N}
    """
    Compute the flux Jacobian for the DVM equations in 1D velocity space.
    """
    return Diagonal(equations.c_vec)
end

function Trixi.max_abs_speed_naive(f_l, f_r, orientation::Integer, equations::DVMEquations1D)
    """
    Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
    """
    λ_l = Trixi.max_abs_speeds(f_l, equations)
    λ_r = Trixi.max_abs_speeds(f_r, equations)
    λ_max = max(λ_l, λ_r)
    return λ_max
end

function Trixi.max_abs_speeds(f, equations::DVMEquations1D)
    """
    Estimate the flux Jacobian eigenvalues by means of Gerschgorin
    """
    return maximum(abs.(real.(eigen(flux_jacobian(f, equations)).values)))
end


# ========================================================================================== #

# Initial Conditions
struct InitialConditionsDVM{N}
    """
    Initial conditions for the DVM equations in 1D velocity space.
    Setups up a Riemann problem with left and right states.
    """
    left::SVector{N}
    right::SVector{N}

    function InitialConditionsDVM(f_left, f_right, equations::DVMEquations1D{N}) where {N}
        left = f_left.(equations.c_vec) #? ones(SVector{N}) * f_left
        right = f_right.(equations.c_vec) #? ones(SVector{N}) * f_right
        @assert all(left .>= 0.0) && all(right .>= 0.0)

        return new{N}(left, right)
    end
end

function (ic::InitialConditionsDVM)(coords, t, equations::DVMEquations1D{N}) where {N}
    if coords[1] < 0.0; return ic.left; else; return ic.right; end
end




# ========================================================================================== #

# physical variables
# u_i = ∫ c^i f(c) dc, i=1,...,N
# u_i(x) = ∫ c_k^i f(c_k)|_x dc, k=1,...,N
function ρ_v_θ_p_DVM(semi, sol, equations::DVMEquations1D)
    """
    Extract physical variables from the DVM solution.

    # Arguments:
    - `semi`: Semi-discretization object (Trixi.jl).
    - `sol`: Solution object containing the solution at different time steps (Trixi.jl).
    - `equations::DVMEquations1D`: DVM equations object (DVMEquations1D).

    # Returns:
    - `x::Vector{Float64}`: Spatial coordinates.
    - `ρ::Vector{Float64}`: Density at each spatial coordinate.
    - `v::Vector{Float64}`: Velocity at each spatial coordinate.
    - `Θ::Vector{Float64}`: Temperature at each spatial coordinate.
    - `p::Vector{Float64}`: Pressure at each spatial coordinate.
    """
    _, coords, variables = collect1DTreeArrays_local(semi, sol.u[end], cons2cons)
    x = coords[2:end-1]
    Fmat = variables[2:end-1, :]'  # First column is density
    ρ = dc(equations) * sum(Fmat, dims=1)
    v = dc(equations) * sum(Fmat .* equations.c_vec, dims=1) ./ ρ
    Θ = dc(equations) * sum(Fmat .* (equations.c_vec .- v).^ 2, dims=1) ./ ρ
    p = ρ .* Θ
    q = dc(equations) * sum(Fmat .* (equations.c_vec .- v).^ 3, dims=1)
    return vec(x), vec(ρ), vec(v), vec(Θ), vec(p), vec(q)
end

function plot_ρ_v_p_DVM(ρ, v, p, x; xlims=(-2.0, 2.0))
    """
    Plot physical variables (density, velocity, pressure) from DVM solution.

    # Arguments:
    - `ρ::Vector{Float64}`: Density array.
    - `v::Vector{Float64}`: Velocity array.
    - `p::Vector{Float64}`: Pressure array.
    - `x::Vector{Float64}`: Spatial coordinates array.
    - `xlims::Tuple{Float64, Float64}=(-2.0, 2.0)`: Tuple specifying x-axis limits for the plot.

    # Returns:
    - `plt`: Plot object containing the plotted variables.
    """
    # create primary plot and plot density and pressure on it
    plt = plot(x, ρ;
        xlabel = "x",
        label = "ρ",
        color = :blue,
        legend = :topleft,
    )

    # plot pressure on the same (primary) axis
    plot!(plt, x, p;
        label = "p",
        color = :green,
    )

    # create a twin y-axis for velocity
    ax2 = twinx(plt)
    plot!(ax2, x, v;
        label = "v",
        color = :red,
        legend = :topright,
    )

    # set x limits on the primary axis (affects both)
    xlims!(plt, xlims)

    return plt
end


# setup for 1D Riemann
function setupDVM1DRiemann(
    N, Kn,
    c_l, c_u, 
    ic_left, ic_right;
    base_tree_level = 8,
    surface_flux = flux_lax_friedrichs,
    volume_flux = flux_central,
    polydeg = 1,        # DG polynomial degree
    domain = (-5.0, 5.0),
    )
    """
    Helper function to set the 1D DVM equations with Riemann initial conditions up.

    # Arguments
    - `N`: Number of discrete velocity points.
    - `Kn`: Knudsen number.
    - `c_l`: Lower bound of the velocity domain.
    - `c_u`: Upper bound of the velocity domain.
    - `ic_left`: Function defining the left initial condition.
    - `ic_right`: Function defining the right initial condition.
    - `base_tree_level=8`: Base refinement level for the mesh.
    - `surface_flux`=flux_lax_friedrichs: Numerical flux function for the surface integral.
    - `volume_flux`=flux_central: Numerical flux function for the volume integral.

    # Returns
    - `basis`: DG basis functions.
    - `mesh`: Computational mesh.
    - `equations`: DVM equations object.
    - `initial_condition`: Initial conditions object.
    - `solver`: DGSEM solver object.
    - `boundary_conditions`: Boundary conditions object.
    """
    equations = DVMEquations1D(N, c_l, c_u, Kn)
    initial_condition = InitialConditionsDVM(
        ic_left,
        ic_right,
        equations
    )

    #= set up semidiscretization =#
    basis = LobattoLegendreBasis(polydeg)

    # ==================================== volume integral treatment ==================================== #
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

    # ==============================================  Solver  ============================================ #
    solver = DGSEM(basis, surface_flux, volume_integral)

    mesh = TreeMesh(
        (domain[1],), 
        (domain[2],), 
        initial_refinement_level=base_tree_level, 
        n_cells_max=10_000, 
        periodicity=false
    )

    boundary_conditions = (
        x_neg = BoundaryConditionDirichlet(initial_condition), 
        x_pos = BoundaryConditionDirichlet(initial_condition)
    )


    return basis, mesh, equations, initial_condition, solver, boundary_conditions
end
