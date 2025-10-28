using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Parameter
N = 50 # ! 250
c_l = -6.0
c_u = 6.0
Kn = 1.0 # Knudsen number
source = relaxation_source

domain = (-5.0, 5.0)
T_end = 0.3

ρ_L = 7.0; v_L = 0.0; θ_L = 1.0
ρ_R = 1.0; v_R = 0.0; θ_R = 1.0

basis, mesh, equations, initial_condition, solver, boundary_conditions = setupBGK1DRiemann(
    N, Kn,
    c_l, c_u, 
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R);
    base_tree_level = 8,
    domain = domain,
)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    boundary_conditions=boundary_conditions, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="Riemann1D/bgk" #! "Riemann1D/bgk_solution_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)", # name of output files
)

#= solve =#
sol = solve(
    ode, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    callback = callbacks,
    saveat = range(tspan[1], tspan[2], length=100)
);

summary_callback()


# Post Processing
# Plotting density, velocity, 
x = LinRange(domain[1], domain[2], length(sol.u[end]) ÷ N)
ρ, v, θ, p = ρ_v_θ_p_BGK(sol.u[end], equations)

p1 = plot_ρ_v_p_bgk(ρ, v, p, x; xlims=(-2.0, 2.0))
display(p1)
savefig(p1, "out/Riemann1D/bgk_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.pdf")

# Store primitive variables in CSV file
CSV.write(
    "out/Riemann1D/bgk_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)


# Analysis of the macroscopic variables compared to the output of the macroscopic variables from the Maxwellian distribution
f_Maxwellian = [Maxwellian(ρ, v, θ).(equations.c_vec) for (ρ, v, θ) in zip(ρ, v, θ)]
ρ_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1]) for i in 1:size(f_Maxwellian, 1)] # density
v_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1] .* equations.c_vec) / ρ_Maxwellian[i] for i in 1:size(f_Maxwellian, 1)] # velocity
θ_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1] .* (equations.c_vec .- v_Maxwellian[i]).^ 2) / ρ_Maxwellian[i] for i in 1:size(f_Maxwellian, 1)] # temperature
p_Maxwellian = ρ_Maxwellian .* θ_Maxwellian # pressure

plt = plot(x, ρ;
    xlabel = "x",
    label = "ρ",
    color = :blue,
    legend = :topleft,
)
plot!(
    plt, x, ρ_Maxwellian;
    label = "ρ (Maxwellian)",
    linestyle = :dash,
    color = :blue,
)
# plot pressure on the same (primary) axis
plot!(plt, x, p;
    label = "p",
    color = :green,
)
plot!(
    plt, x, p_Maxwellian;
    label = "p (Maxwellian)",
    linestyle = :dash,
    color = :green,
)


# create a twin y-axis for velocity
ax2 = twinx(plt)
plot!(ax2, x, v;
    label = "v",
    color = :red,
    legend = :topright,
)
plot!(
    ax2, x, v_Maxwellian;
    label = "v (Maxwellian)",
    linestyle = :dash,
    color = :red,
)

display(plt)

# Error-Analysis (of the final solution)
L2_error_ρ = norm(ρ .- ρ_Maxwellian, 2) / norm(ρ_Maxwellian, 2)
L2_error_v = norm(v .- v_Maxwellian, 2) / norm(v_Maxwellian, 2)
L2_error_p = norm(p .- p_Maxwellian, 2) / norm(p_Maxwellian, 2)

println("L2 relative error in density: $L2_error_ρ")
println("L2 relative error in velocity: $L2_error_v")
println("L2 relative error in pressure: $L2_error_p")

# ========================================================================================== #
# Error-Analysis over the time-steps
L2_error_ρ = []; L2_error_v = []; L2_error_p = [];
for i in 1:size(sol.u, 1)
    ρ, v, θ, p = ρ_v_θ_p_BGK(sol.u[i], equations)

    f_Maxwellian = [Maxwellian(ρ, v, θ).(equations.c_vec) for (ρ, v, θ) in zip(ρ, v, θ)]
    ρ_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1]) for i in 1:size(f_Maxwellian, 1)] # density
    v_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1] .* equations.c_vec) / ρ_Maxwellian[i] for i in 1:size(f_Maxwellian, 1)] # velocity
    θ_Maxwellian = [HyQMOM.dc(equations) * sum(f_Maxwellian[i,:][1] .* (equations.c_vec .- v_Maxwellian[i]).^ 2) / ρ_Maxwellian[i] for i in 1:size(f_Maxwellian, 1)] # temperature
    p_Maxwellian = ρ_Maxwellian .* θ_Maxwellian # pressure

    push!(L2_error_ρ, norm(ρ .- ρ_Maxwellian, 2) / norm(ρ_Maxwellian, 2))
    push!(L2_error_v, norm(v .- v_Maxwellian, 2) / norm(v_Maxwellian, 2))
    push!(L2_error_p, norm(p .- p_Maxwellian, 2) / norm(p_Maxwellian, 2))
end

plt = plot(
    sol.t[3:end], L2_error_ρ[3:end];
    xlabel = "t",
    label = "ρ",
    color = :blue,
    legend = :topleft,
)
plot!(
    plt, sol.t[3:end], L2_error_v[3:end];
    label = "v",
    color = :red,
)
plot!(
    plt, sol.t[3:end], L2_error_p[3:end];
    label = "p",
    color = :green,
)
plot!(plt, yaxis=:log)
display(plt)
