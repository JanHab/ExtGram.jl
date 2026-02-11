using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LaTeXStrings

# Parameter
N = 50
c_l = -6.0
c_u = 6.0
Kn = 1.0 # ! doesn't matter, as zero-relaxation in the vlasov_poisson_source_term # Knudsen number
T_end = 15.0
source = vlasov_poisson_source
x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

base_tree_level = 5

# todo: need to export BGKEquations1D
equations = ExtGram.BGKEquations1D(N, c_l, c_u, Kn)
ρ0 = 1.0; v0 = 0.0; θ0 = 1.0
ϵ = 0.001
k = 0.5
# todo: need to export InitialConditionsLandauDamping_BGK
initial_condition = ExtGram.InitialConditionsLandauDamping_BGK(
    ρ0,
    ϵ,
    v0,
    θ0,
    k,
    equations
)

#= set up semidiscretization =#
polydeg = 1 #!1 #!1
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

# boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition), x_pos = BoundaryConditionDirichlet(initial_condition))

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    # boundary_conditions=boundary_conditions, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

# callbacks, summary_callback = callbacksGramianMomentEquations(
#     semi, tspan, basis; 
#     cfl = 0.99,          # Maximum cfl number
#     plot_interval = 20,  # plot every 20 steps
#     time_interval = 500, # save at 500 time intervals
#     # * name="VlasovPoisson/VlasovPoisson_moments_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level)", # name used for output
#     name="VlasovPoisson/gram_moments", # name used for output
# )
# # Add Vlasov-Poisson callback
# callbacks = CallbackSet(callbacks, vlasov_poisson_callback(;M, mesh, domain))
# todo: need to export vlasov_poisson_callback_BGK
callbacks = CallbackSet(ExtGram.vlasov_poisson_callback_BGK(;N=N, mesh, domain, equations))

#= solve =#
sol = solve(
    ode, 
    # Euler();
    CarpenterKennedy2N54(
        williamson_condition = false
    );
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., 
    save_everystep=true,
    callback = callbacks,
);

summary_callback()

x = LinRange(domain[1], domain[2], length(sol.u[end]) ÷ N)
ρ, v, θ, p, q = ρ_v_θ_p_BGK(sol.u[end], equations)

p1 = plot_ρ_v_p_bgk(ρ, v, p, x; xlims=(domain[1], domain[2]))
display(p1)
# savefig(p1, "out/Riemann1D/bgk_Kn$(Kn)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_ρ_v_p.pdf")

# Access the energy history
E_L2_history = ExtGram.ELECTRIC_FIELD_BGK.E_L2
time_callback = ExtGram.ELECTRIC_FIELD_BGK.times  # Use actual times from callback

# You can then plot it or analyze it
γ = -0.1533 # theoretical decay rate for k=1/2
γt = exp.(γ .* time_callback)
plot(
    time_callback, E_L2_history ./ E_L2_history[1],
    xlabel="Time", 
    label="BGK N=$N (from callback)",
    ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
    yaxis=:log
)
plot!(
    time_callback, γt,
    label="Theoretical Decay exp($γ t)", 
    linestyle=:dash
)
plot!(legend=:bottomleft)
savefig("out/VlasovPoisson/BGK_energy_from_callback.pdf")
# store to csv file
CSV.write(
    "out/VlasovPoisson/BGK_energy_moments_T$(T_end)_N$(N)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level).csv",
    Tables.columntable((
        time=time_callback, E_L2=E_L2_history, 
        E_L2_normalized=E_L2_history ./ E_L2_history[1], 
        theoretical_decay=γt
    ))
)
