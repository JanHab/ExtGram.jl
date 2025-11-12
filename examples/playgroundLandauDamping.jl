using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
using LaTeXStrings

# Parameter
M = 4    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0 # ! doesn't matter, as zero-relaxation in the vlasov_poisson_source_term # Knudsen number
T_end = 15.0
source = vlasov_poisson_source
x_lower = 0.0; x_upper = 4.0*π
domain = (x_lower, x_upper)

base_tree_level = 8 #!5 # ! 8

equations = GramianMomentEquations1D(M, Kn, closure)
ρ0 = 1.0; v0 = 0.0; θ0 = 1.0
ϵ = 0.001
k = 0.5
initial_condition = InitialConditionsLandauDamping(
    ρ0,
    ϵ,
    v0,
    θ0,
    k,
    equations
)

#= set up semidiscretization =#
polydeg = 3 #!1 #!1
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


plot(HyQMOM.ELECTRIC_FIELD.x_range[2:end-1], HyQMOM.ELECTRIC_FIELD.E[2:end-1];
    xlabel="x", ylabel="E(x)", label="Electric Field at t=$T_end", lw=2,
)
plot(HyQMOM.ELECTRIC_FIELD.x_range[2:end-1], HyQMOM.ELECTRIC_FIELD.ρ[2:end-1];
    label="ρ", lw=2, linestyle=:dash    
)

u = sol.u[end]
MP1 = HyQMOM.ELECTRIC_FIELD.MP1
n_cells = length(u) ÷ MP1

# Reshape to extract density
U_matrix = reshape(u, MP1, n_cells)
ρ = U_matrix[1, :]  # First row is density

plot(ρ;
    label="ρ (from solution)", lw=2, linestyle=:dot, marker=:x
)
plot!(
    HyQMOM.ELECTRIC_FIELD.variables[2:end-1,1], marker=:o
)

# Access the energy history
E_L2_history = HyQMOM.ELECTRIC_FIELD.E_L2
time_callback = HyQMOM.ELECTRIC_FIELD.times  # Use actual times from callback

# You can then plot it or analyze it
γ = -0.1533 # theoretical decay rate for k=1/2
γt = exp.(γ .* time_callback)
plot(
    time_callback, E_L2_history ./ E_L2_history[1],
    xlabel="Time", 
    label="HyQMOM M=$M (from callback)",
    ylabel=L"∥E(t,⋅)∥_{L^2} / ∥E(0,⋅)∥_{L^2}", 
    yaxis=:log
)
plot!(
    time_callback, γt,
    label="Theoretical Decay exp($γ t)", 
    linestyle=:dash
)
plot!(legend=:bottomleft)

# nodes = basis.nodes
# weights = basis.weights

# n_elements = length(mesh.tree)
# x = Matrix{Float64}(undef, length(nodes), n_elements)
# coordinates_min = domain[1]
# coordinates_max = domain[2]
# dx = (coordinates_max - coordinates_min) / n_elements
# for element in 1:n_elements
#     x_l = coordinates_min + (element - 1) * dx + dx / 2
#     for i in eachindex(nodes)
#         ξ = nodes[i] # nodes in [-1, 1]
#         x[i, element] = x_l + dx / 2 * ξ
#     end
# end

# println(size(E))
# println(size(x))
# x_ = vec(x)
# # unique_x = vec(unique(x_))
# unique_x = vec(unique(x_ -> round(x_, digits=4), x))
# println(size(unique_x))
# s = Int(size(ode.u0, 1)/(M+1))
# sol_end = sol.u[end]
# U_matrix_end = reshape(sol_end, M+1, length(sol_end) ÷ (M+1))
# ρ_end = U_matrix_end[1, :]  # First row is density

# x_vector = []
# if size(x, 1) == 1
#     for i in 1:size(x, 1)
#         push!(x_vector, x[i, 1])
#     end
# else
#     for i in 1:size(x, 2)
#         if i == 1
#             for j in 1:size(x, 1)
#                 push!(x_vector, x[j, i])
#             end
#         # elseif i == size(x, 2)
#         #     for j in 1:size(x, 1)
#         #         push!(x_vector, x[j, i])
#         #     end
#         else
#             for j in 2:size(x, 1)
#                 push!(x_vector, x[j, i])
#             end
#         end
#     end
# end
# x_vector

# U_matrix = reshape(ode.u0, M+1, length(ode.u0) ÷ (M+1))
# ρ = U_matrix[1, :]  # First row is density

# plot(vec(x), x -> 3 * x^2, label = "f'", lw = 2)
# scatter!(vec(x), x -> 3 * x^2, label = "f'", lw = 2)
