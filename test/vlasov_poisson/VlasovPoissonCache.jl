# Tests for the Vlasov-Poisson workspace: the polydeg guard, the FFT Poisson solve and the field lookup.
using Trixi, Test

L = 4π
mesh = TreeMesh((0.0,), (L,), initial_refinement_level = 3, n_cells_max = 100, periodicity = true)
ncells = 2^3

@testset "polydeg 1 and 2 give an equispaced vertex grid" begin
    for p in (1, 2)
        solver = DGSEM(LobattoLegendreBasis(p), flux_lax_friedrichs)
        cache = VlasovPoissonCache(mesh, solver)
        @test length(cache.x) == p * ncells + 1
        h = cache.x[2] - cache.x[1]
        @test all(d -> isapprox(d, h; rtol = 1e-12), diff(cache.x))
        @test cache.x[end] - cache.x[1] ≈ L
        @test cache.x₁ == cache.x[1]
        @test cache.h ≈ h
        @test cache.n == p * ncells == length(cache.ρ)
    end
end

@testset "polydeg ≥ 3 is rejected" begin
    for p in (3, 4, 5)
        solver = DGSEM(LobattoLegendreBasis(p), flux_lax_friedrichs)
        @test_throws ArgumentError VlasovPoissonCache(mesh, solver)
    end
end

@testset "Poisson solve reproduces E = (ε/k) sin(kx) for ρ = 1 + ε cos(kx)" begin
    ε, k = 0.01, 0.5
    for p in (1, 2)
        solver = DGSEM(LobattoLegendreBasis(p), flux_lax_friedrichs)
        cache = VlasovPoissonCache(mesh, solver)
        cache.ρ .= 1 .+ ε .* cos.(k .* cache.x[1:end-1])
        ExtGram.solve_poisson!(cache)
        E_exact = ε / k .* sin.(k .* cache.x)
        @test maximum(abs.(cache.E .- E_exact)) < 1e-12 * maximum(abs.(E_exact))
        @test cache.E[end] == cache.E[1]
        # field lookup: exact at every vertex, periodic image at the last one, linear in between
        @test all(j -> ExtGram.field_at(cache, cache.x[j]) ≈ cache.E[j], eachindex(cache.x))
        @test ExtGram.field_at(cache, cache.x[end]) == cache.E[1]
        x_mid = 0.5 * (cache.x[3] + cache.x[4])
        @test ExtGram.field_at(cache, x_mid) ≈ 0.5 * (cache.E[3] + cache.E[4])
    end
end
