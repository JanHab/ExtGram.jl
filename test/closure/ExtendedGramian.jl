closure = "ExtGram"
Kn = 0.5    # irrelevant

# todo: update to new closure
# @testset "M3_extended" begin
#     M = 3

#     equations = GramianMomentEquations1D(M, Kn, closure)
    
#     u = SVector{M+1,Float64}(rand(M+1)...).*10

#     trixi_closure = HyQMOM.closure(u, equations)
#     correct_closure = test_closure(u, equations)
#     @test trixi_closure ≈ correct_closure
# end

@testset "M4_extended" begin
    M = 4

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end

# todo: update to new closure
# @testset "M5_extended" begin
#     M = 5

#     equations = GramianMomentEquations1D(M, Kn, closure)
    
#     u = SVector{M+1,Float64}(rand(M+1)...).*10

#     trixi_closure = HyQMOM.closure(u, equations)
#     correct_closure = test_closure(u, equations)
#     @test trixi_closure ≈ correct_closure
# end

@testset "M6_extended" begin
    M = 6

    equations = GramianMomentEquations1D(M, Kn, closure)
    
    u = SVector{M+1,Float64}(rand(M+1)...).*10

    trixi_closure = HyQMOM.closure(u, equations)
    correct_closure = test_closure(u, equations)
    @test trixi_closure ≈ correct_closure
end

@testset "Compare with Mathematica for M=3" begin
    testmoments_1 = [1.1, 2.2, 3.1, 4.4]
    testmoments_2 = [3.5, 100.3, 22.1, -10.2]
    testmoments_3 = [12_313.123123, 433.123213, 22.34, 10_000]
    testmoments = [testmoments_1, testmoments_2, testmoments_3]
    M = 3

    # Classical
    equations = GramianMomentEquations1D(M, Kn, "Gram")
    mathematica_closure = [6.24406, -5.64747, 1.40732e+7]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end

    # todo: update to new closure
    # Extended
    # equations = GramianMomentEquations1D(M, Kn, "ExtGram")
    # mathematica_closure = [6.93077, -401.198, 1.05553e+7]
    # for i in eachindex(testmoments)
    #     trixi_closure = closure(testmoments[i], equations)
    #     @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    # end
end

@testset "Compare with Mathematica for M=4" begin
    testmoments_1 = [1.1, 2.2, 3.1, 4.4, 6.4]
    testmoments_2 = [3.5, 100.3, 22.1, -10.2, 22.4]
    testmoments_3 = [12313.123123, 433.123213, 22.34, 1000, 11.213]
    testmoments = [testmoments_1, testmoments_2, testmoments_3]
    M = 4

    # Classical
    equations = GramianMomentEquations1D(M, Kn, "Gram")
    mathematica_closure = [9.07692, 6.59831, -3368.4]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end

    # Extended
    equations = GramianMomentEquations1D(M, Kn, "ExtGram")
    mathematica_closure = [9.40081, 16.0905, -2.96487e+7]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end
end

@testset "Compare with Mathematica for M=5" begin
    testmoments_1 = [1.1, 2.2, 3.1, 4.4, 7.2, 3.4]
    testmoments_2 = [3.5, 100.3, 22.1, -10.2, 102.2, 234.3]
    testmoments_3 = [12_313.123123, 433.123213, 22.34, 10_000, 2342.2, 239.3]
    testmoments = [testmoments_1, testmoments_2, testmoments_3]
    M = 5

    # Classical
    equations = GramianMomentEquations1D(M, Kn, "Gram")
    mathematica_closure = [63.7832, 382.892, 8116.65]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end

    # todo: update to new closure
    # Extended
    # equations = GramianMomentEquations1D(M, Kn, "ExtGram")
    # mathematica_closure = [37.1617, 309.211, -150923.]
    # for i in eachindex(testmoments)
    #     trixi_closure = closure(testmoments[i], equations)
    #     @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    # end
end

@testset "Compare with Mathematica for M=6" begin
    testmoments_1 = [1.1, 2.2, 3.1, 4.4, 7.2, 3.4, 7.3]
    testmoments_2 = [3.5, 100.3, 22.1, -10.2, 102.2, 234.3, 723.2]
    testmoments_3 = [12_313.123123, 433.123213, 22.34, 10_000, 2342.2, 239.3, 234.32]
    testmoments = [testmoments_1, testmoments_2, testmoments_3]
    M = 6

    # Classical
    equations = GramianMomentEquations1D(M, Kn, "Gram")
    mathematica_closure = [-0.635406, 1408.32, 1944.15]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end

    # Extended
    equations = GramianMomentEquations1D(M, Kn, "ExtGram")
    mathematica_closure = [533.868, 2290.58, -148.104]
    for i in eachindex(testmoments)
        trixi_closure = HyQMOM.closure(testmoments[i], equations)
        @test isapprox(trixi_closure, mathematica_closure[i], atol=1e-4, rtol=1e-4)
    end
end
