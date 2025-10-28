if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, StaticArrays

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
