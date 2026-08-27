
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays

M = 4
Kn = 1.0
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
base_tree_level = 6 # ! Increase later
source = vlasov_maxwell_source 
T_end = 20.0
cfl = 0.99 # ? Decrease later? 0.45
time_interval = 10 # ? Increase later? 10
name = "out"

Bx, By, Bz = 0.0, 0.0, 0.0 # * Vanishing magnetic contributions for the moment
v1, v2, v3, α, k = 0.0, 0.0, 0.0, 0.5, 0.4
domain = (0.0, 2.0*π/k)

equations = GramianMomentEquations1D3V(M, Kn, "ExtGram")

initial_condition = InitialConditionsVlasovMaxwellLandauDamping(
    v1, v2, v3, α, k,
    equations
)

basis = LobattoLegendreBasis(polydeg)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=true)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=cfl)

callbacks = CallbackSet(
    alive_callback, summary_callback, stepsize_callback,
    vlasov_maxwell_callback(;Bx, By, Bz, M, mesh, domain)
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
    # saveat = range(tspan[1], tspan[2], length=100)
);

plot(ExtGram.VLASOV_MAXWELL_FIELD.times, ExtGram.VLASOV_MAXWELL_FIELD.Ex_L2, xlabel="Time", ylabel="L2-norm of E_x", title="L2-norm of Electric Field over Time")
savefig("L2_norm_Electric_Field.png")

#= ---------- initial state: ρ, v, p, θ over x ---------- =#
u_end = Array(Trixi.wrap_array(sol.u[1], semi))          # (nvars, nnodes, nelements)
x     = vec(Array(semi.cache.elements.node_coordinates))   # same (node, element) ordering

nvars = nvariables(equations)
prim  = vec([cons2prim(SVector{nvars}(u_end[:, i, e]), equations)
             for i in axes(u_end, 2), e in axes(u_end, 3)])

i200, i020, i002 = equations._pressure_index      # (5, 8, 10) for M = 4

ρ   = [p[1] for p in prim]
v_x   = [p[2] for p in prim]
v_y = [p[3] for p in prim]
v_z = [p[4] for p in prim]
θ   = [(p[i200] + p[i020] + p[i002]) / (3 * p[1]) for p in prim]
prs = ρ .* θ                                       # scalar pressure p = ρθ

perm = sortperm(x); xs = x[perm]
plt = plot(
    plot(xs, ρ[perm], ylabel="ρ", title="density"),
    plot(xs, [v_x[perm] v_y[perm] v_z[perm]], label=["v_x" "v_y" "v_z"], title="velocity", legend=:topright),
    plot(xs, prs[perm], ylabel="p", title="pressure"),
    plot(xs, θ[perm], ylabel="θ", title="temperature"),
    layout=(2,2), size=(900,600), xlabel="x", lw=2, legend=false,
    plot_title="1D-3V  M=$M  Kn=$Kn  t=$(round(sol.t[end], digits=3))"
)
display(plt)
savefig("Initial_State.png")


#= ---------- Final state: ρ, v, p, θ over x ---------- =#
u_end = Array(Trixi.wrap_array(sol.u[end], semi))          # (nvars, nnodes, nelements)
x     = vec(Array(semi.cache.elements.node_coordinates))   # same (node, element) ordering

nvars = nvariables(equations)
prim  = vec([cons2prim(SVector{nvars}(u_end[:, i, e]), equations)
             for i in axes(u_end, 2), e in axes(u_end, 3)])

i200, i020, i002 = equations._pressure_index      # (5, 8, 10) for M = 4

ρ   = [p[1] for p in prim]
v_x   = [p[2] for p in prim]
v_y = [p[3] for p in prim]
v_z = [p[4] for p in prim]
θ   = [(p[i200] + p[i020] + p[i002]) / (3 * p[1]) for p in prim]
prs = ρ .* θ                                       # scalar pressure p = ρθ

perm = sortperm(x); xs = x[perm]
plt = plot(
    plot(xs, ρ[perm], ylabel="ρ", title="density"),
    plot(xs, [v_x[perm] v_y[perm] v_z[perm]], label=["v_x" "v_y" "v_z"], title="velocity", legend=:topright),
    plot(xs, prs[perm], ylabel="p", title="pressure"),
    plot(xs, θ[perm], ylabel="θ", title="temperature"),
    layout=(2,2), size=(900,600), xlabel="x", lw=2, legend=false,
    plot_title="1D-3V  M=$M  Kn=$Kn  t=$(round(sol.t[end], digits=3))"
)
display(plt)
savefig("Final_State.png")