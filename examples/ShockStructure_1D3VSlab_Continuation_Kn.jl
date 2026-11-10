
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays

M = parse(Int, ARGS[1])
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = parse(Float64, ARGS[3])
source_string = ARGS[4]
source = source_string == "relaxation_source" ? relaxation_source : zero_source 
angle_name = ARGS[5] # "Max", "Arc", "Det", "Bad"
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
domain = (-30.0, 30.0)
base_tree_level = 10
T_end = 50.0
cfl = 0.45
time_interval = 10 # ? Increase later? 10
slab_geometry = true # ? false


angles = get(SlabAngles, "M$(M)_$(angle_name)", nothing)
    if angles === nothing
        error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", \"Con\".")
    end
theta = [angles[i][1] for i in 1:length(angles)]
phi = [angles[i][2] for i in 1:length(angles)]

Ma = parse(Float64, ARGS[6]) # 3.8 and 6.5 in Cai paper
rho_L = 1.0
rho_R = (4*Ma^2) / (Ma^2 + 3)
v1_L, v2_L, v3_L = sqrt(5/3) * Ma, 0, 0
v1_R, v2_R, v3_R = sqrt(5/3) * (Ma^2 + 3) / (4 * Ma), 0, 0
theta_L = 1.0
theta_R = (5*Ma^2 - 1) * (Ma^2 + 3) / (16 * Ma^2)
f_left = Maxwellian1D3V(rho_L, (v1_L, v2_L, v3_L), theta_L)
f_right = Maxwellian1D3V(rho_R, (v1_R, v2_R, v3_R), theta_R)

equations = GramianMomentEquations1D3V(M, 0.01, closure, slab_geometry=slab_geometry, theta=theta, phi=phi)

initial_condition = InitialConditionsShockTube1D3V(
    f_left, # Density, velocity, temperature
    f_right, # Shock in density, but not velocity, temperature initially
    M,
    equations
)

basis = LobattoLegendreBasis(polydeg)

# shock capturing
indicator_sc = IndicatorHennemannGassner(
    equations, basis,
    variable = (u, eqns)->u[2]^2/u[1] # * Placeholder for the moment
)

volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((domain[1],), (domain[2],), initial_refinement_level=base_tree_level, n_cells_max=10_000, periodicity=false)

boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition), x_pos = BoundaryConditionDirichlet(initial_condition))

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    boundary_conditions=boundary_conditions, 
    source_terms=source
)

t_end_vec = [10.0, 10.0, 10.0, 30.0] # ToDo: Adapt the numbers
Kn_vec = [0.01, 0.1, 0.5, 1.0]
for i in 1:length(t_end_vec)
    Kn = Kn_vec[i]
    T_end = t_end_vec[i]

    equations = GramianMomentEquations1D3V(M, Kn, closure, slab_geometry=slab_geometry, theta=theta, phi=phi)

    semi.equations = equations

    #= set up ODE =#
    tspan = (0.0, T_end)
    ode = semidiscretize(semi, tspan)
    if i > 1
        ode.u0 .= u_prev
    end

    alive_callback = AliveCallback(analysis_interval=100)
    summary_callback = SummaryCallback()
    stepsize_callback = StepsizeCallback(cfl=cfl)

    name = "1D3V/1D3V_ShockStructure/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_anglepair$(angle_name)_MaL$(Ma)"
    save_solution_cons = SaveTriangulationCallback(
        time_interval=tspan[2]/time_interval,
        save_initial_solution=true,
        file_format="tsv",
        append_solution=true,
        solution_variables = cons2cons,
        clear_out_dir=false,
        name=name * "_cons",
        info="basis = $(Base.typename(typeof(basis)).wrapper)"
    )
    save_solution_prim = SaveTriangulationCallback(
        time_interval=tspan[2]/time_interval,
        save_initial_solution=true,
        file_format="tsv",
        append_solution=true,
        solution_variables = cons2prim,
        clear_out_dir=false,
        name=name * "_prim",
        info="basis = $(Base.typename(typeof(basis)).wrapper)"
    )
    callbacks = CallbackSet(
        alive_callback, summary_callback, stepsize_callback,
        save_solution_cons, save_solution_prim
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
    u_prev = copy(sol.u[end])
end



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
    plot_title="1D-3V  M=$M  Kn=$Kn  t=$(round(sol.t[end], digits=3))",
    xlims=(-10, 10)
)
display(plt)
# savefig("Final_State.png")

sigma_xx = [p[i200] for p in prim] # σ_xx
q_x = [p[i200] for p in prim] # Q_xx
plt = plot(
    plot(xs, sigma_xx[perm], ylabel="σ_xx", title="σ_xx"),
    plot(xs, q_x[perm], ylabel="Q_xx", title="Q_xx"),
    layout=(1,2), size=(900,400), xlabel="x", lw=2, legend=false,
    plot_title="1D-3V  M=$M  Kn=$Kn  t=$(round(sol.t[end], digits=3))",
    xlims=(-10,10)
)