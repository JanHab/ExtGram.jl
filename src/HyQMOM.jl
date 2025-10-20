module HyQMOM

# External dependencies
using LinearAlgebra, StaticArrays, Trixi, FastGaussQuadrature, ForwardDiff, OrdinaryDiffEq, Plots, Interpolations, Statistics

# Gramian Moments Implementation
include("gramian_moment_equations.jl")
include("testing.jl")
include("initial_conditions.jl")

export GramianMomentEquations1D, gramian, closure, moment_prim2cons, moment_cons2prim, dCdu, flux_jacobian, relaxation_source, zero_source, numerical_flux
export InitialConditionsShockTube, Maxwellian, convective_moments, primitive_moments
export test_closure, check_realizability

# Vlasov-Poisson implementation
# implemented via a callback
# The Electric fields’ contribution is added in the source term. The electric field is calculated after each time step in a Discrete Callback with a Jacobi iteration.
# ϕ^i = 1/2 (ρ^i (Δx)^2 + ϕ^{i-1} + ϕ^{i+1})
# The boundary conditions are periodic.
# This means we get: ϕ^0 = ϕ^N and ϕ^{N+1} = ϕ^1
# In our case: ϕ^0 = ϕ^{N+1} = 0
# The electric field is then calculated via E = -∂ϕ/∂x ≈ -(ϕ^{i+1} - ϕ^{i-1})/(2Δx)
# needed dependencies: StaticArrays, Interpolations
using FFTW
include("vlasov_poisson.jl")
export vlasov_poisson_callback, vlasov_poisson_source, InitialConditionsCosine

# TrixiTree2Triangulation
# This module implements the triangulation of a 2D TreeMesh with subsequent
# writing of the mesh and solution into a basic file readable by external visualization programs.
# The intended use is as callback and is therefore simple to use.
# Required packages: Trixi, OrdinaryDiffEq, HDF5, WriteVTK, (LinearAlgebra, Printf)
include("TrixiTree2Triangulation/TrixiTree2Triangulation.jl")
export SaveTriangulationCallback

# Semidiscretization and callbacks
include("setup1D.jl")
export setupGramianMomentEquations1DRiemann, callbacksGramianMomentEquations

# analysis tools
include("analysis.jl")
export plot_ρ_v_p, plot_λ_max, TVD_space, readsol, readfile

# BGK equation
include("bgk_equation.jl")
export setupBGK1DRiemann, InitialConditionsBGK, ρ_v_θ_p_BGK, setupBGK1DRiemann, plot_ρ_v_p_bgk

end # module HyQMOM
