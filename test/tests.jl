if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, StaticArrays

using Test

@testset "GradClosure" begin
    include("closure/Grad.jl")
end

@testset "GramianClosure" begin
    include("closure/Gramian.jl")
end

@testset "ExtendedGramianClosure" begin
    include("closure/ExtendedGramian.jl")
end

@testset "1D3 - Matrix Transformations" begin
    include("1D3V/TestTransformations.jl")
end

@testset "1D3V - full" begin
    include("1D3V/1D3V_full.jl")
end

@testset "1D3V - Reduction to Slab-geometry" begin
    include("1D3V/1D3V_slab.jl")
end