using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

######################################################
################## Velocity Variation ################
######################################################
M_vector = [5, 9, 13, 17, 25, 51]
Kn = 1.0 # Doesn't matter for the closure, necessary for defining the equations
closures = ["Gram", "ExtGram", "Grad"]

ρ1, v1, θ1 = 0.4, 0.0, 1.0
ρ2, θ2 = 0.6, 1.0

v2_vector = range(0.5, 4.0, length=50)

colors = [:blue, :red, :green]
n = length(M_vector)
pl = Vector{Any}(undef, n)# = scatter(layout=4)
for i in 1:n
    pl[i] = scatter()
end
for (i, M) in enumerate(M_vector)
    for closure in closures
        L2 = Float64[]
        for v2 in v2_vector
            equations = GramianMomentEquations1D(M, Kn, closure)

            f1 = Maxwellian(ρ1, v1, θ1)
            f2 = Maxwellian(ρ2, v2, θ2)
            momentList = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))

            nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
            push!(L2, abs((nextMoment - momentList[end]) / momentList[end]))
        end

        scatter!(
            pl[i],
            v2_vector, L2, 
            yscale=:log10, 
            label="",
            # label=closure, 
            title="M = $M",
            xlabel="v2", ylabel="ϵᵣ", 
            color=colors[findfirst(==(closure), closures)],
        )
    end
end
# labels = ["Gram", "ExtGram", "Grad"]
labels = ["\n  $(closures[1])\n" "\n  $(closures[2])\n" "\n  $(closures[3])\n"]
n_labels = length(closures)
# colors = palette(:default)[1:n_labels]'
p0 = scatter(
    (-n_labels:-1)',(-n_labels:-1)', 
    lims=(0,1), legendfontsize=7, 
    legend=:left, fg_color_legend = nothing, 
    label=labels, fc=colors, 
    frame=:none
);
# display(pl)
l = @layout  [grid(3,2) a{0.2w}]
plot(pl..., p0, layout=l, size=(800,500))
savefig("out/MorinClosure/ErrorVsVelocity.pdf")

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
        L2 = Float64[]
        for M in M_vector
            equations = GramianMomentEquations1D(M, Kn, closure)

            f1 = Maxwellian(ρ1, v1, θ1)
            f2 = Maxwellian(ρ2, v2, θ2)
            momentList = convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))

            nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
            push!(L2, abs((nextMoment - momentList[end]) / momentList[end]))
        end

        plot!(
            pl[i],
            M_vector, L2, 
            yscale=:log10, 
            label="",
            marker=:o,
            title="v₂ = $v2",
            xlabel="M", ylabel="ϵᵣ", 
            color=colors[findfirst(==(closure), closures)],
            # legend=:bottomright;
        )
    end
end
# display(pl)
l = @layout  [grid(2,2) a{0.2w}]
plot(pl..., p0, layout=l, size=(800,500))
savefig("out/MorinClosure/ErrorVsM.pdf")


######################################################
################## Equilibrium Test ##################
######################################################
M_vector = 3:2:101
ρ, v, θ = 1.0, 0.0, 1.0

pl = plot()
for (i, closure) in enumerate(closures)
    L2 = Float64[]
    for M in M_vector
        equations = GramianMomentEquations1D(M, Kn, closure)
        f = Maxwellian(ρ, v, θ)
        momentList = convective_moments(f, Val(M+2))
        nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
        push!(L2, abs((nextMoment - momentList[end]) / momentList[end]))
    end

    plot!(
        pl,
        M_vector, L2, 
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
savefig("out/MorinClosure/EquilibriumPreservation.pdf")
