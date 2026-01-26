using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables

using Trapz

f0(v, ϕ, v0) = 1 / √(2π) * exp(-.5 * (sign(v) * √(v^2 - 2ϕ) - v0)^2)
f1(v, ϕ, v0, β) = 1 / √(2π) * exp(-.5 * (β  * (v^2 - 2ϕ) + v0)^2)
f(v, ϕ, v0, β) = if v^2 > 2ϕ; f0(v, ϕ, v0); else; f1(v, ϕ, v0, β); end;


# Parameter
ψ = 0.2
Δ = 1
x0 = 0.0
ϕ(x) = ψ * exp(((x-x0)/Δ)^2) #10.
v0 = 1.
β = -0.1

x = LinRange(-4, 4, 500)
v = LinRange(-5, 5, 100)

heatmap(x, v, [f(vj, ϕ(xi), v0, β) for xi in x, vj in v]', c=:bluesreds, xlabel="x", ylabel="v", title="Electron hole distribution function f(x,v)", colorbar_title="f(x,v)")


x = 0.0
plot(v, f.(v, ϕ(x), v0, β), xlabel="v", ylabel="f(v)", title="Electron hole distribution function at x=0", legend=false)


# Parameter
M = 4    # number of moments
closure = "Gram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.3
source = zero_source #!relaxation_source

# Domain and discretization parameters
x_lower = -4.0; x_upper = 4.0 #! -4.0, 4.0
domain = (x_lower, x_upper)
polydeg = 1  # polynomial degree
base_tree_level = 9 #!8  # initial mesh refinement level
surface_flux = flux_lax_friedrichs
volume_flux = flux_central

equations = GramianMomentEquations1D(M, Kn, closure)

# Parameter
ψ = 0.1 # 0.2
Δ = 1#!4
v0 = 0.1 #!1.
β = -0.1
initial_condition = InitialConditionsElectronHole(
    ElectronHole(ψ, Δ, v0, β), # Density, velocity, temperature
    equations
)

#= set up semidiscretization =#
basis = LobattoLegendreBasis(polydeg)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=false) # ! do we use a periodic domain here?

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

######### Plot initial condition moments #########
# Plot initial distribution function
x_coords = LinRange(x_lower, x_upper, 500)
# v_coords = LinRange(-5.0, 5.0, 100)
# f_vals = [initial_condition((x, ), 0.0, equations) for x in x_coords, v in v_coords]
# heatmap(x_coords, v_coords, [f_vals[i,j][1] for i in 1:length(x_coords), j in 1:length(v_coords)]', xlabel="x", ylabel="v", title="Electron hole initial distribution function f(x,v)", colorbar_title="f(x,v)")

# Plot for different x-values
# x_slices = [-2.0, 0.0, 2.0]
# plt = plot()
# for x in x_slices
#     plot!(plt, v_coords, [initial_condition((x, ), 0.0, equations)[1] for v in v_coords], label="x = $x", xlabel="v", ylabel="f(v)", title="Electron hole initial distribution function slices at different x")
# end
# display(plt)

# Plot initial moments
_, coords, vars = ExtGram.collect1DTreeArrays_local(semi, ode.u0, cons2cons)
x = coords[2:end-1]  # exclude ghost cells
ρ = vars[2:end-1, 1]
v = vars[2:end-1, 2] ./ ρ
θ = (vars[2:end-1, 3] ./ ρ) .- v.^2

plt = plot()
plot!(plt, x, ρ, label="ρ (initial)", xlabel="x", ylabel="Density / Velocity / Temperature", title="Electron hole initial moments")
plot!(plt, x, v, label="v (initial)")
plot!(plt, x, θ, label="θ (initial)")
display(plt)

# = Callbacks =#
callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.9,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
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
x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)
display(p1)