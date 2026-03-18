#=
Author: Matthias Geratz | Matthias.Geratz@rwth-aachen.de
March 2023

This module implements an OrdinaryDiffEq callback for the triangulation of a 2D TreeMesh with subsequent
writing of the mesh and solution into a basic file readable by external visualization programs.
Note that the current procedure degrades solution quality quite a bit, due to triangulating only the element corner nodes,
solution averaging at shared nodes and linear output basis functions.
This is a proof of concept code and therefore may exhibit bugs and bad performance. User caution is advised.

Required packages: Trixi, OrdinaryDiffEq, HDF5, WriteVTK, (LinearAlgebra, Printf)

Installation from this local directory:
using Pkg
Pkg.develop(PackageSpec(path = "<path to module folder>/TrixiTree2Triangulation"))
=#


# module TrixiTree2Triangulation # commented out to make ExtGram installation easier

using OrdinaryDiffEq: DiscreteCallback, u_modified!
using DiffEqCallbacks: PeriodicCallback, PeriodicCallbackAffect
using Trixi, LinearAlgebra
using Printf, HDF5, WriteVTK


# 2D
include("converter.jl")
include("writer.jl")
# 1D
include("converter_1d.jl")
include("writer_1d.jl")

include("methods.jl")
include("callback.jl") # callback definition happens here

# export SaveTriangulationCallback

# end # module TrixiTree2Triangulation
