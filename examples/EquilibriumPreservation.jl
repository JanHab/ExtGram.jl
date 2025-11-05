using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# M_vector = [3, 5, 7, 9, 11, 13]    # number of moments
M_vector = Int[]
for M=3:2:101
    push!(M_vector, M)
end
closure = "ExtGram"
Kn = 1.0 # Doesn't matter for the closure, necessary for defining the equations

ρ, v, θ = 1.0, 0.0, 1.0

for M in M_vector
    equations = GramianMomentEquations1D(M, Kn, closure)

    f = Maxwellian(ρ, v, θ)
    momentList = convective_moments(f, Val(M+2))

    nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
    if !isapprox(nextMoment, momentList[end]; rtol=1e-6)
        println("Testing closure for M = $M moments")
        println("Next moment predicted by closure: $nextMoment")
        println("Exact next moment: $(momentList[end])")
        println("Relative error: $((nextMoment - momentList[end]) / momentList[end])")
        # error("Closure did not reproduce the correct next moment for M = $M")
    end

end