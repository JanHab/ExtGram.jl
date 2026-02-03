using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

# Parameter
M = 8    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.3 #!25.0 #!0.3
source = relaxation_source

# Domain and discretization parameters
x_lower = -2.0; x_upper = 2.0
# x_lower = -20.0; x_upper = 100.0
# x_lower = -50; x_upper = 50
domain = (x_lower, x_upper)
polydeg = 1  # polynomial degree
base_tree_level = 8 #!8  # initial mesh refinement level
surface_flux = flux_lax_friedrichs
volume_flux = flux_central

ρ_L = 7.0
v_L = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R = 0.0
θ_R = 1.0

# Rankine-Hugoniot conditions for shock tube
# ρ_L = 1.0
# θ_L = 1.0
# p_L = ρ_L * θ_L
# Ma = 1.0001
# γ = 5.0/3.0
# ρ_R = ρ_L * (Ma^2 * (γ+1)) / (2 + Ma^2 * (γ - 1))
# p_R = p_L * (1 - γ + 2 * γ * Ma^2) / (1 + γ)
# θ_R = p_R / ρ_R
# c_L = sqrt(γ * p_L / ρ_L) # speed of sound left
# v_L = Ma * c_L
# v_R = v_L * ρ_L / ρ_R

# Cai
# Ma = 1.4

# ρ_L = 1.0
# ρ_R = 2*Ma^2 / (Ma^2 + 1)

# v_L = sqrt(3) * Ma
# v_R = sqrt(3) / 2 * (Ma^2 + 1) / Ma

# θ_L = 1.0
# θ_R = (3*Ma^2 - 1) * (Ma^2 + 1) / (4 * Ma^2)


# Setting up everything
equations = GramianMomentEquations1D(M, Kn, closure)
initial_condition = InitialConditionsShockTube(
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R), # Shock in density, but not velocity, temperature initially
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)

# shock capturing
# #= crashes for p > 1, i.e. this does not help at all
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    alpha_max = 0.5, #  α_max = 1.0 seems natural -> corresponds to pure first order FV (Gassner paper)
    alpha_min = 0.001,
    alpha_smooth = true, #* false, # smoothes with all neighboring indicators to remove numerical artifacts
    variable = (u, eqns)->u[1]*u[3]
) # ? Seems to restrict to M>=4
# `custom variable for smoothness detection?`

volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)
# volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh(
    (domain[1],), (domain[2],), 
    initial_refinement_level=base_tree_level, 
    n_cells_max=10_000, 
    periodicity=false
)

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
    cfl = 0.9,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
)

# amr_controller = ControllerThreeLevel(
#     semi, indicator_sc;
#     base_level=base_tree_level-3,
#     med_level=base_tree_level, med_threshold=0.1,
#     max_level=base_tree_level+3, max_threshold=0.6
# )

# amr_callback = AMRCallback(
#     semi, amr_controller,
#     interval=5,
#     adapt_initial_condition=true,
#     adapt_initial_condition_only_refine=true
# )

# callbacks = CallbackSet(callbacks, amr_callback)



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
x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)

θ = p ./ ρ

plot(x, ρ[:, end], label="ρ", lw=2)
plot(x, v[:, end], label="v", lw=2)
plot(x, p[:, end], label="p", lw=2)
plot(x, θ[:, end], label="θ", lw=2)

# Eigenvalues
u_final = sol.u[end]
L = length(u_final)
Nloc = L ÷ (M+1)
Fmat = reshape(u_final, M+1, Nloc)

x_index_center = Nloc ÷ 2
u0 = Fmat[1, x_index_center]
u1 = Fmat[2, x_index_center]
u2 = Fmat[3, x_index_center]
u3 = Fmat[4, x_index_center]
u4 = Fmat[5, x_index_center]
u5 = Fmat[6, x_index_center]
u6 = Fmat[7, x_index_center]
u7 = Fmat[8, x_index_center]
u8 = Fmat[9, x_index_center]


u = SVector{M+1}(u0, u1, u2, u3, u4, u5, u6, u7, u8)
jacobian = flux_jacobian(u, equations)
using LinearAlgebra
λ = real.(eigen(jacobian).values)
y = zeros(length(λ))

# plot values of eigenvalues on a horizontal line
plot(λ, y, seriestype=:scatter, title="Eigenvalues at center cell", xlabel="Index", ylabel="Eigenvalue")

# Plot primitive variables at final time
# plot(x, (ρ[:, end] .- ρ_L) ./ (ρ_R - ρ_L), label="ρ", lw=2, xlim=(-10, 10))
# plot(x, (v[:, end] .- v_R) ./ (v_L - v_R), label="v", lw=2, xlim=(-10, 10))
# plot(x, (p[:, end] .- θ_L .* ρ_L) ./ (θ_R .* ρ_R - θ_L .* ρ_L), label="p", lw=2, xlim=(-10, 10))
# plot(x, (θ[:, end] .- θ_L) ./ (θ_R - θ_L), label="θ", lw=2, xlim=(-10, 10))

pl = plot(xlim=(-10, 10), title="T=$(T_end)", size=(500,500), yticks=0:0.1:1);
plot!(pl, x, (ρ[:, end] .- ρ_L) ./ (ρ_R - ρ_L), label="ρ", lw=2);
plot!(pl, x, (v[:, end] .- v_R) ./ (v_L - v_R), label="v", lw=2);
# plot!(pl, x, (p[:, end] .- θ_L .* ρ_L) ./ (θ_R .* ρ_R - θ_L .* ρ_L), label="p", lw=2)
plot!(pl, x, (θ[:, end] .- θ_L) ./ (θ_R - θ_L), label="θ", lw=2);
display(pl)


display(p1)
savefig(p1, "out/Riemann1D/ρ_v_p.pdf")

# Store primitive variables in CSV file
CSV.write(
    "out/Riemann1D/ρ_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)

# Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
n_plots = 5
p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
# display(p2)
savefig(p2, "out/Riemann1D/λ_max.pdf")

# Plot total variation in space over time
p3, mass, momentum, energy = conservation(sol, M, semi)
# display(p3)
savefig(p3, "out/Riemann1D/conservation.pdf")
