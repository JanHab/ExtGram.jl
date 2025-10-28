closure = "Gram"
Kn = 0.5    # irrelevant

@testset "M3_standard" begin
    M = 3

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end

@testset "M4_standard" begin
    M = 4

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end

@testset "M5_standard" begin
    M = 5

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end

@testset "M6_standard" begin
    M = 6

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end