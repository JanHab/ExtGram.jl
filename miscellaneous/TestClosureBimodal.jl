using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

######################################################
################## Velocity Variation ################
######################################################
M_vector = [4, 5, 8, 9, 12, 13, 16, 17, 24, 25, 50, 51]
Kn = 1.0 # Doesn't matter for the closure, necessary for defining the equations
closures = ["Gram", "ExtGram", "Grad"]

ρ1, v1, θ1 = 0.4, 0.0, 1.0
ρ2, θ2 = 0.6, 1.0

v2_vector = range(0.5, 4.0, length=200)

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

            nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
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
# labels = ["Gram", "ExtGram", "Grad"]
# labels = ["\n  $(closures[1])\n" "\n  $(closures[2])\n" "\n  $(closures[3])\n"]
# n_labels = length(closures)
# # colors = palette(:default)[1:n_labels]'
# p0 = scatter(
#     (-n_labels:-1)',(-n_labels:-1)', 
#     lims=(0,1), legendfontsize=7, 
#     legend=:left, fg_color_legend = nothing, 
#     label=labels, fc=colors, 
#     frame=:none
# );
# # display(pl)
# l = @layout  [grid(3,2) a{0.2w}]
# plot(pl..., p0, layout=l, size=(800,500))
# savefig("out/Figures/CompareClosure/ErrorVsVelocity.pdf")

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

            nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
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
# display(pl)
# l = @layout  [grid(2,2) a{0.2w}]
# plot(pl..., p0, layout=l, size=(800,500))
# savefig("out/Figures/CompareClosure/ErrorVsM.pdf")

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

            nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
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
# display(pl)
# l = @layout  [grid(2,2) a{0.2w}]
# plot(pl..., p0, layout=l, size=(800,500))
# savefig("out/Figures/CompareClosure/ErrorVsM.pdf")

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
        nextMoment = HyQMOM.closure(momentList[1:end-1], equations)
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

















using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end 
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, QuadGK

######################################################
################## Distribution Helpers ##############
######################################################

# 1. Mott-Smith Moments
function get_mott_smith_moments(M, x, Ma, gamma=5/3)
    # Rankine-Hugoniot downstream values
    rho_star = (Ma^2 * (gamma + 1)) / (2 + Ma^2 * (gamma - 1))
    v_star = Ma / rho_star
    theta_star = (1 - gamma + 2*gamma*Ma^2) / ((1 + gamma) * rho_star)
    v_scale = sqrt(gamma)
    
    z = 1 / (1 + exp(x))
    
    # Components as Maxwellians
    f1 = Maxwellian(z, Ma * v_scale, 1.0)
    f2 = Maxwellian((1-z)*rho_star, v_star * v_scale, theta_star)
    
    # Return moments up to order M+1
    return convective_moments(f1, Val(M+2)) + convective_moments(f2, Val(M+2))
end

# 2. Electron-Hole Moments (requires numerical integration for the trapped part)
function get_eh_moments(M, phi, v0, beta)
    moments = zeros(M+2)
    const_factor = (2*pi)^(-0.5)
    
    for k in 0:(M+1)
        f_eh = v -> begin
            v2 = v^2
            threshold = 2*phi
            if v2 > threshold
                # Untrapped
                return const_factor * exp(-0.5 * (sign(v)*sqrt(v2 - threshold) - v0)^2)
            else
                # Trapped
                return const_factor * exp(-0.5 * (beta*(v2 - threshold) + v0^2))
            end
        end
        val, _ = quadgk(v -> (v^k) * f_eh(v), -Inf, Inf, rtol=1e-10)
        moments[k+1] = val
    end
    return moments
end

######################################################
################## Global Parameters #################
######################################################
M_vector = [3, 4, 5, 8] # Small sample for quick testing
closures = ["Gram", "ExtGram", "Grad"]
colors = [:blue, :red, :green]
Kn = 1.0

# Legend Placeholder
# p_legend = scatter(
#     (1:3)', (1:3)'#, mc=colors'#, label=permutedims(closures), 
#     # frame=:none, grid=false, showaxis=false
# )

######################################################
################## Mott-Smith Analysis ###############
######################################################
println("Running Mott-Smith Analysis...")
Ma = 3.0
x_range = range(-5.0, 5.0, length=100)
pl_ms = [plot(title="M=$M", xlabel="x", ylabel="ϵᵣ", yscale=:log10) for M in M_vector]

for (i, M) in enumerate(M_vector)
    for (c_idx, closure) in enumerate(closures)
        eqs = GramianMomentEquations1D(M, Kn, closure)
        errs = Float64[]
        for x in x_range
            m_list = get_mott_smith_moments(M, x, Ma)
            next_m = HyQMOM.closure(m_list[1:end-1], eqs)
            push!(errs, abs((next_m - m_list[end]) / m_list[end]))
        end
        plot!(pl_ms[i], x_range, errs, color=colors[c_idx], label="")
    end
end
l_ms = @layout [grid(2,2) a{0.15w}]
plot(pl_ms..., layout=l_ms, size=(900, 500)) |> display
# savefig("Accuracy_MottSmith_Ma$(Ma).pdf")

######################################################
################## Electron-Hole Analysis ############
######################################################
println("Running Electron-Hole Analysis...")
phi_range = range(0.1, 2.0, length=50)
v0_eh, beta_eh = 1.5, -0.05
pl_eh = [plot(title="M=$M", xlabel="ϕ", ylabel="ϵᵣ", yscale=:log10) for M in M_vector]

for (i, M) in enumerate(M_vector)
    for (c_idx, closure) in enumerate(closures)
        eqs = GramianMomentEquations1D(M, Kn, closure)
        errs = Float64[]
        for phi in phi_range
            m_list = get_eh_moments(M, phi, v0_eh, beta_eh)
            next_m = HyQMOM.closure(m_list[1:end-1], eqs)
            push!(errs, abs((next_m - m_list[end]) / m_list[end]))
        end
        plot!(pl_eh[i], phi_range, errs, color=colors[c_idx], label="")
    end
end
l_eh = @layout [grid(2,2) a{0.15w}]
plot(pl_eh..., layout=l_eh, size=(900, 500)) |> display
# savefig("Accuracy_ElectronHole_beta$(beta_eh).pdf")
