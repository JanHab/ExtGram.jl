# reference values are taken from a predefined run of the Trixi example to make sure that changes to the closure do not affect results unexpectedly
# The values were compared to a Mathematica implementation of the Grad closure beforehand and shown only slight deviations due to different numerical implementations

atol = 1e-12
rtol = 1e-12

Kn = 0.5 # doesn't matter for closure test
closure = "Grad"

ρ1 = 0.4
v1 = 0.0
θ1 = 1.0
ρ2 = 0.6
v2 = [0.0, 0.5, 1.0]
θ2 = 1.0

f1 = Maxwellian(ρ1, v1, θ1)
f2_vec = [Maxwellian(ρ2, v_2, θ2) for v_2 in v2]

function trixi_closure(f1, f2, M, Kn, closure)
    u = convective_moments(f1, Val(M+1)) .+ convective_moments(f2, Val(M+1))
    equations = GramianMomentEquations1D(M, Kn, closure)

    return HyQMOM.closure(u, equations) # todo: remove HyQMOM prefix
end

@testset "M4" begin
    M = 4
    @test isapprox(trixi_closure(f1, f2_vec[1], M, Kn, closure), 0.0, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[2], M, Kn, closure), 5.265929999999987, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[3], M, Kn, closure), 15.50975999999995, atol=1e-12, rtol=1e-12)
end

@testset "M5" begin
    M = 5
    @test isapprox(trixi_closure(f1, f2_vec[1], M, Kn, closure), 15.000000000000231, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[2], M, Kn, closure), 22.318845000000323, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[3], M, Kn, closure), 51.406080000000756, atol=1e-12, rtol=1e-12)
end

@testset "M6" begin
    M = 6
    @test isapprox(trixi_closure(f1, f2_vec[1], M, Kn, closure), 0.0, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[2], M, Kn, closure), 39.774802500000376, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[3], M, Kn, closure), 139.3747200000013, atol=1e-12, rtol=1e-12)
end

@testset "M7" begin
    M = 7
    @test isapprox(trixi_closure(f1, f2_vec[1], M, Kn, closure), 105., atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[2], M, Kn, closure), 176.14127756999986, atol=1e-12, rtol=1e-12)
    @test isapprox(trixi_closure(f1, f2_vec[3], M, Kn, closure), 500.76705792000155, atol=1e-12, rtol=1e-12)
end
