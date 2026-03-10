if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LaTeXStrings

# Access arguments by index
M = parse(Int, ARGS[1])
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
T_end = parse(Float64, ARGS[3])
base_tree_level = parse(Int, ARGS[4]) # e.g. 8
polydeg = parse(Int, ARGS[5]) # e.g. 3
ρ_0 = parse(Float64, ARGS[6]) # 1.0
v_0 = parse(Float64, ARGS[7]) # 0.0
θ_0 = parse(Float64, ARGS[8]) # 1.0
ϵ = parse(Float64, ARGS[9]) # 0.01
k = parse(Float64, ARGS[10]) # 0.5

x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

Kn = 1.0 # doesn't matter, as zero-relaxation in the vlasov_poisson_source_term
source = vlasov_poisson_source


equations = GramianMomentEquations1D(M, Kn, closure)
initial_condition = InitialConditionsLandauDamping(
    ρ_0,
    ϵ,
    v_0,
    θ_0,
    k,
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

surface_flux = flux_lax_friedrichs
volume_flux = flux_central
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=true)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    # boundary_conditions=boundary_conditions, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.99,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    time_interval = 500, # save at 500 time intervals
    name="VlasovPoisson/LandauDamping/moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name used for output
)
# Add Vlasov-Poisson callback
callbacks = CallbackSet(callbacks, vlasov_poisson_callback(;M, mesh, domain))

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
E_L2_history = ExtGram.ELECTRIC_FIELD.E_L2
time_callback = ExtGram.ELECTRIC_FIELD.times  # Use actual times from callback

# You can then plot it or analyze it
# plot(
#     time_callback, E_L2_history ./ E_L2_history[1],
#     xlabel="Time", 
#     label="ExtGram M=$M (from callback)",
#     ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
#     yaxis=:log
# )
# if k == 0.5
#     γ = -0.1533 # theoretical decay rate for k=1/2
#     γt = exp.(γ .* time_callback)
#     plot!(
#         time_callback, γt,
#         label="Theoretical Decay exp($γ t)", 
#         linestyle=:dash
#     )
# end
# plot!(legend=:bottomleft)
# savefig("out/VlasovPoisson/LandauDamping/energy_moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_p$(polydeg)_level$(base_tree_level)_x_lower$(x_lower)_x_upper$(x_upper).pdf")
# store to csv file
CSV.write(
    "out/VlasovPoisson/LandauDamping/energy_moments_solution_M$(M)_closure$(closure)_Kn$(Kn)_T_end$(T_end)_rho_0$(ρ_0)_v_0$(v_0)_theta_0$(θ_0)_epsilon$(ϵ)_k$(k)_p$(polydeg)_level$(base_tree_level)_x_lower$(x_lower)_x_upper$(x_upper).csv",
    Tables.columntable((
        time=time_callback, E_L2=E_L2_history, 
        E_L2_normalized=E_L2_history ./ E_L2_history[1], 
    ))
)