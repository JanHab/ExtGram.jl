using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LaTeXStrings

# Parameter
M = 6    # number of moments
closure = "Grad" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0 # ! doesn't matter, as zero-relaxation in the vlasov_poisson_source_term # Knudsen number
T_end = 15.0
source = vlasov_poisson_source
x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

base_tree_level = 5

equations = GramianMomentEquations1D(M, Kn, closure)
ρ0 = 1.0; v0 = 0.0; θ0 = 1.0
ϵ = 0.001
k = 0.5
initial_condition = InitialConditionsCosine(
    ρ0,
    ϵ,
    v0,
    θ0,
    k,
    equations
)

#= set up semidiscretization =#
polydeg = 3
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
    # * name="VlasovPoisson/VlasovPoisson_moments_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level)", # name used for output
    name="VlasovPoisson/gram_moments", # name used for output
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
E_L2_history = HyQMOM.ELECTRIC_FIELD.E_L2;
time_callback = HyQMOM.ELECTRIC_FIELD.times;  # Use actual times from callback

# You can then plot it or analyze it
γ = -0.1533; # theoretical decay rate for k=1/2
γt = exp.(γ .* time_callback);
plot(
    time_callback, E_L2_history ./ E_L2_history[1],
    xlabel="Time", 
    label="HyQMOM M=$M (from callback)",
    ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
    yaxis=:log,
    legend=:bottomleft
)
plot!(
    time_callback, γt,
    label="Theoretical Decay exp($γ t)", 
    linestyle=:dash
)
savefig("out/VlasovPoisson/energy_from_callback.pdf")
# store to csv file
CSV.write(
    "out/VlasovPoisson/energy_moments_T$(T_end)_M$(M)_k$(k)_ϵ$(ϵ)_p$(polydeg)_level$(base_tree_level).csv",
    Tables.columntable((
        time=time_callback, E_L2=E_L2_history, 
        E_L2_normalized=E_L2_history ./ E_L2_history[1], 
        theoretical_decay=γt
    ))
)
