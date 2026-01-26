using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = 8
Kn = 1.0 # Doesn't matter here
closure = "ExtGram"

# Define two Maxwellians
ρ1, v1, θ1 = 0.4, 0.0, 1.0
ρ2, v2, θ2 = 0.6, 2.0, 1.0

# Equations
equations = GramianMomentEquations1D(M, Kn, closure)

# Moments
f1 = Maxwellian(ρ1, v1, θ1)
f2 = Maxwellian(ρ2, v2, θ2)

# conservative
momentListConservative = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))
nextMomentConservative = ExtGram.closure(momentListConservative[1:end-1], equations)

# primitive
momentListPrimitive = ExtGram.moment_cons2prim(momentListConservative, equations)
nextMomentPrimitive = ExtGram.closure(momentListPrimitive[1:end-1], equations)
primtive_to_conservative = ExtGram.moment_prim2cons(SVector{M+2,Float64}(momentListPrimitive[1:end-1]..., nextMomentPrimitive), equations)

# Check conversations
println("moments(conservative): $momentListConservative")
println("moments(primitive -> conservative): $primtive_to_conservative")

# Check accuracy
println("Closure(conservative): $nextMomentConservative")
println("Exact(conservative): $(momentListConservative[end])")
println("ϵ_r(conservative): $(abs((nextMomentConservative - momentListConservative[end]) / momentListConservative[end]))")
println("ϵ_r(primitive -> conservative): $(abs((primtive_to_conservative[end] - momentListConservative[end]) / momentListConservative[end]))")
println("Closure(primitive): $nextMomentPrimitive")
println("Exact(primitive): $(momentListPrimitive[end])")
println("ϵ_r(primitive): $(abs((nextMomentPrimitive - momentListPrimitive[end]) / momentListPrimitive[end]))")



v2_vector = range(0.5, 4.0, length=200)
ϵ_rel_conservative = Float64[]
ϵ_rel_transformation = Float64[]

for v2 in v2_vector
    f2 = Maxwellian(ρ2, v2, θ2)

    # conservative
    momentListConservative = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))
    nextMomentConservative = ExtGram.closure(momentListConservative[1:end-1], equations)
    push!(ϵ_rel_conservative, abs((nextMomentConservative - momentListConservative[end]) / momentListConservative[end]))

    # primitive
    momentListPrimitive = ExtGram.moment_cons2prim(momentListConservative, equations)
    nextMomentPrimitive = ExtGram.closure(momentListPrimitive[1:end-1], equations)
    primtive_to_conservative = ExtGram.moment_prim2cons(SVector{M+2,Float64}(momentListPrimitive[1:end-1]..., nextMomentPrimitive), equations)
    nextMomentTransformation = primtive_to_conservative[end]
    push!(ϵ_rel_transformation, abs((nextMomentTransformation - momentListConservative[end]) / momentListConservative[end]))
end

scatter(
    v2_vector, 
    [ϵ_rel_conservative ϵ_rel_transformation], 
    yscale=:log10, 
    label=["conservative" "closure(primitive) -> conservative"], 
    title="M = $M, closure = $closure",
    xlabel="v₂", ylabel="ϵᵣ", 
)