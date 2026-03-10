using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

######################################################
################## Velocity Variation ################
######################################################
M_vector = [4, 5, 8, 9, 12, 13, 16, 17, 24, 25, 50, 51]
Kn = 1.0 # Doesn't matter for the closure, necessary for defining the equations
closures = ["Gram", "ExtGram", "Grad"]

ρ1, v1, θ1 = 0.4, 0.0, 1.0
ρ2, θ2 = 0.6, 1.0

v2_vector = range(1e-4, 4.0, length=200)

colors = [:blue, :red, :green]
n = length(M_vector)
pl = Vector{Any}(undef, n)# = scatter(layout=4)
for i in 1:n
    pl[i] = scatter()
end
for (i, M) in enumerate(M_vector)
    for closure in closures
        ϵ_rel = Float64[]
        for v2 in v2_vector
            equations = GramianMomentEquations1D(M, Kn, closure)

            f1 = Maxwellian(ρ1, v1, θ1)
            f2 = Maxwellian(ρ2, v2, θ2)
            momentList = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))

            nextMoment = ExtGram.closure(momentList[1:end-1], equations)
            push!(ϵ_rel, abs((nextMoment - momentList[end]) / momentList[end]))
        end

        scatter!(
            pl[i],
            v2_vector, ϵ_rel, 
            yscale=:log10, 
            label="",
            # label=closure, 
            title="M = $M",
            xlabel="v2", ylabel="ϵᵣ", 
            color=colors[findfirst(==(closure), closures)],
        )

        CSV.write(
            "out/CompareClosure/Accuracy_over_v2_closure_$(closure)_M$(M).csv",
            Tables.columntable((
                v2_vector=v2_vector, rel_err=ϵ_rel
            ))
        )
    end
end

######################################################
################## Order of Moments M ################
######################################################
v2_vector = [2.0, 3.0, 4.0, 5.0]
M_vector = 5:2:25

n = length(v2_vector)
pl = Vector{Any}(undef, n)# = scatter(layout=4)
for i in 1:n
    pl[i] = plot()
end
for (i, v2) in enumerate(v2_vector)
    for closure in closures
        ϵ_rel = Float64[]
        for M in M_vector
            equations = GramianMomentEquations1D(M, Kn, closure)

            f1 = Maxwellian(ρ1, v1, θ1)
            f2 = Maxwellian(ρ2, v2, θ2)
            momentList = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))

            nextMoment = ExtGram.closure(momentList[1:end-1], equations)
            push!(ϵ_rel, abs((nextMoment - momentList[end]) / momentList[end]))
        end

        plot!(
            pl[i],
            M_vector, ϵ_rel, 
            yscale=:log10, 
            label="",
            marker=:o,
            title="v₂ = $v2",
            xlabel="M", ylabel="ϵᵣ", 
            color=colors[findfirst(==(closure), closures)],
            # legend=:bottomright;
        )

        CSV.write(
            "out/CompareClosure/Accuracy_over_OddMoments_closure_$(closure)_v2$(v2).csv",
            Tables.columntable((
                M_vector=M_vector, rel_err=ϵ_rel
            ))
        )
    end
end

M_vector = 4:2:24

n = length(v2_vector)
pl = Vector{Any}(undef, n)# = scatter(layout=4)
for i in 1:n
    pl[i] = plot()
end
for (i, v2) in enumerate(v2_vector)
    for closure in closures
        ϵ_rel = Float64[]
        for M in M_vector
            equations = GramianMomentEquations1D(M, Kn, closure)

            f1 = Maxwellian(ρ1, v1, θ1)
            f2 = Maxwellian(ρ2, v2, θ2)
            momentList = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))

            nextMoment = ExtGram.closure(momentList[1:end-1], equations)
            push!(ϵ_rel, abs((nextMoment - momentList[end]) / momentList[end]))
        end

        plot!(
            pl[i],
            M_vector, ϵ_rel, 
            yscale=:log10, 
            label="",
            marker=:o,
            title="v₂ = $v2",
            xlabel="M", ylabel="ϵᵣ", 
            color=colors[findfirst(==(closure), closures)],
            # legend=:bottomright;
        )

        CSV.write(
            "out/CompareClosure/Accuracy_over_EvenMoments_closure_$(closure)_v2$(v2).csv",
            Tables.columntable((
                M_vector=M_vector, rel_err=ϵ_rel
            ))
        )
    end
end

######################################################
################## Equilibrium Test ##################
######################################################
M_vector = 3:2:101
ρ, v, θ = 1.0, 0.0, 1.0

pl = plot()
for (i, closure) in enumerate(closures)
    ϵ_rel = Float64[]
    for M in M_vector
        equations = GramianMomentEquations1D(M, Kn, closure)
        f = Maxwellian(ρ, v, θ)
        momentList = convective_moments(f, Val(M+2))
        nextMoment = ExtGram.closure(momentList[1:end-1], equations)
        push!(ϵ_rel, abs((nextMoment - momentList[end]) / momentList[end]))
    end

    plot!(
        pl,
        M_vector, ϵ_rel, 
        yscale=:log10, 
        label=closure,
        marker=:o,
        title="Equilibrium Preservation",
        xlabel="M", ylabel="ϵᵣ", 
        color=colors[i],
    )
end

plot!(
    M_vector, ones(length(M_vector)).*1e-6, 
    yscale=:log10, 
    label="1e-6",
    linestyle=:dash,
    color=:black,
)

plot!(
    pl, legend=:left
)

display(pl)
# savefig("out/CompareClosure/EquilibriumPreservation.pdf")
