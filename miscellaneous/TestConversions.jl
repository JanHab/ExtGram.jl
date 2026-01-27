using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, Plots

M = 8
closure = "ExtGram"
Kn = 1.0
source = relaxation_source

# IC
ρ1 = 0.4
v1 = 0.0
θ1 = 1.0
# ρ2 = 0.6
ρ2_vector = range(0.1, 4.0, length=200)
v2 = 0.0
# v2_vector = range(0.5, 4.0, length=200)
θ2 = 1.0

equations = GramianMomentEquations1D(M, Kn, closure)

f1 = Maxwellian(ρ1, v1, θ1); # Density, velocity, temperature

ϵ_cons = Float64[]
ϵ_prim = Float64[]
# for v2 in v2_vector
for ρ2 in ρ2_vector
# for v2_shift in v2_vector
    f2 = Maxwellian(ρ2, v2, θ2); # Shock in density, but not velocity, temperature initially
    conservative_moments_full = convective_moments(f1, Val(M+2)) .+ convective_moments(f2, Val(M+2))
    # conservative_moments_full = conservative_moments_full .+ SVector{M+2}(0.0, v2_shift, zeros(M)...)
    conservative_moments = conservative_moments_full[1:end-1] # Remove hidden truth
    conservative_closure = ExtGram.closure(conservative_moments, equations)

    # Print accuracy
    # println("Conservative closure accuracy: ")
    # println("Closure value: ", conservative_closure)
    # println("Hidden truth: ", conservative_moments_full[end])
    # println("Relative error: ", abs(conservative_closure - conservative_moments_full[end]) / abs(conservative_moments_full[end]))
    push!(ϵ_cons, abs((conservative_closure - conservative_moments_full[end]) / conservative_moments_full[end]))


    primitive_moments_full = ExtGram.moment_cons2prim(conservative_moments_full, equations)
    primitive_moments = primitive_moments_full[1:end-1] # Remove hidden truth
    primitive_closure = ExtGram.closure(primitive_moments, equations)
    primitive_moments_full_with_closure = SVector{M+2}(primitive_moments..., primitive_closure)
    conservative_moments_full_with_closure = ExtGram.moment_prim2cons(primitive_moments_full_with_closure, equations)

    # Print accuracy
    # println("Primitive closure accuracy: ")
    # println("Closure value: ", conservative_moments_full_with_closure[end])
    # println("Hidden truth: ", conservative_moments_full[end])
    # println("Relative error: ", abs(conservative_moments_full_with_closure[end] - conservative_moments_full[end]) / abs(conservative_moments_full[end]))
    push!(ϵ_prim, abs((conservative_moments_full_with_closure[end] - conservative_moments_full[end]) / conservative_moments_full[end]))
end

plt = scatter(
    xlabel="ρ2", ylabel="ϵᵣ", yscale=:log10,
    title="M=$M, Closure=$(closure)",
);
scatter!(plt, ρ2_vector, ϵ_cons, label="Conservative", color=:blue);
scatter!(plt, ρ2_vector, ϵ_prim, label="Primitive", color=:red);
display(plt);


f = Vector{Float64}[]
c_range = range(-5.0, 5.0, length=250)
ρ2_range = range(0.1, 4.0, length=5)
for ρ2 in ρ2_range
    f1 = Maxwellian(ρ1, v1, θ1); # Density, velocity, temperature
    f2 = Maxwellian(ρ2, v2, θ2);
    push!(f, f1.(c_range) + f2.(c_range))
end
plt2 = plot(
    # c_range, f[1],
    xlabel="c", ylabel="f(c)",
    title="Distribution functions for different ρ2",
);
for i in 1:length(f)
    plot!(plt2, c_range, f[i], label="ρ2=$(round(ρ2_range[i], digits=2))");
end
display(plt2);