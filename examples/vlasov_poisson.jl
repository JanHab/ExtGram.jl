using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LaTeXStrings

# Parameter
M = 8    # number of moments
extended = true # flag for extended gramian closure or "standard" closure
Kn = 1.0  # Knudsen number
T_end = 25.0
source = vlasov_poisson_source
x_lower = 0.0; x_upper = 2.0*π
domain = (x_lower, x_upper)

base_tree_level = 7

equations = GramianMomentEquations1D(M, Kn, extended)
initial_condition = InitialConditionsCosine(
    1.0, # ρ0
    0.001, # ϵ
    0.0, # v0
    1.0, # θ0
    0.5, # k
    equations
)

#= set up semidiscretization =#
polydeg = 1
basis = LobattoLegendreBasis(polydeg)

# shock capturing
# #= crashes for p > 1, i.e. this does not help at all
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    alpha_max = 1.0,
    alpha_min = 0.01,
    alpha_smooth = true,
    variable = (u, eqns)->u[1]*u[3]
)

surface_flux = flux_lax_friedrichs
volume_flux = flux_central
volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=true) # ! periodic

boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition), x_pos = BoundaryConditionDirichlet(initial_condition))

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
    cfl = 0.99,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="VlasovPoisson_moments",
)
# Add Vlasov-Poisson callback
callbacks = CallbackSet(callbacks, vlasov_poisson_callback(;M, domain))

#= solve =#
sol = solve(
    ode, 
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    save_everystep=true,
    callback = callbacks,
);

summary_callback()

# Access the energy history
E_L2_history = HyQMOM.ELECTRIC_FIELD.E_L2

# You can then plot it or analyze it
# todo: fix to cfl (not fixed time-step)
time = LinRange(0, T_end, length(E_L2_history))
γ = -0.1533 # theoretical decay rate for k=1/2
γt = exp.(γ .* time)
plot(
    time, E_L2_history ./ E_L2_history[1],
    xlabel="Time", 
    label="HyQMOM M=$M",
    ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
    yaxis=:log
)
plot!(
    time, γt,
    label="Theoretical Decay exp($γ t)", 
    linestyle=:dash
)
savefig("out/VlasovPoisson/energy_from_callback.pdf")


L2_history = []
u_final = sol.u[end]
L = length(u_final)
n_cells = L ÷ (M+1)
for j in 1:length(sol.u)
    Fmat = reshape(sol.u[j], M+1, n_cells)
    ρ = Fmat[1, :]
    E = HyQMOM.solve_poisson_periodic_fft(ρ, domain)
    E_L2 = (sum(E.^2) * (domain[2] - domain[1]) / n_cells)^(1/2)
    push!(L2_history, E_L2)
end

plot(
    sol.t, L2_history / L2_history[1], 
    xlabel="Time", ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
    yaxis=:log, 
    label="L2-Norm from Postprocessing"
)
plot!(sol.t, exp.(γ .* sol.t), 
    linestyle=:dash, 
    label="Theoretical Decay exp($γ t)"
)
savefig("out/VlasovPoisson/energy_from_postprocessing.pdf")
