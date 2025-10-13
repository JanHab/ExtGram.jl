module HyQMOM

# External dependencies
using LinearAlgebra, StaticArrays, Trixi, FastGaussQuadrature, ForwardDiff, OrdinaryDiffEq, Plots

# Gramian Moments Implementation
include("gramian_moment_equations.jl")
include("testing.jl")
include("initial_conditions.jl")

export GramianMomentEquations1D, gramian, closure, moment_prim2cons, moment_cons2prim, dCdu, flux_jacobian, relaxation_source, zero_source, numerical_flux
export InitialConditionsShockTube, Maxwellian, convective_moments, primitive_moments
export test_closure, check_realizability

# TrixiTree2Triangulation
# This module implements the triangulation of a 2D TreeMesh with subsequent
# writing of the mesh and solution into a basic file readable by external visualization programs.
# The intended use is as callback and is therefore simple to use.

# Required packages: Trixi, OrdinaryDiffEq, HDF5, WriteVTK, (LinearAlgebra, Printf)

include("TrixiTree2Triangulation/TrixiTree2Triangulation.jl")
export SaveTriangulationCallback

# Semidiscretization and callbacks
include("setup1D.jl")
export setupGramianMomentEquations, callbacksGramianMomentEquations

# analysis tools
include("analysis.jl")
export plot_ρ_v_p, plot_λ_max, TVD_space

end # module HyQMOM
