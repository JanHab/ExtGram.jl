#=
This module implements the triangulation of a 2D TreeMesh with subsequent
writing of the mesh and solution into a basic file readable by external visualization programs.
The intended use is as callback and is therefore simple to use.

Required packages: Trixi, OrdinaryDiffEq, HDF5, WriteVTK, (LinearAlgebra, Printf)

Installation from this local directory:
using Pkg
Pkg.develop(PackageSpec(path = "<path to module folder>/TrixiTree2Triangulation"))
=#


module TrixiTree2Triangulation

include("methods.jl")
include("callback.jl")

export SaveTriangulationCallback

end # module TrixiTree2Triangulation
