
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
angle_name = ARGS[5] # "Max", "Arc", "Det"
polydeg = 1
volume_flux = flux_central
surface_flux = flux_lax_friedrichs
domain = (-2.0, 2.0)
base_tree_level = 5 #!10
T_end = 0.3 
cfl = 0.99 # ? Decrease later? 0.45
time_interval = 10 # ? Increase later? 10
slab_geometry = true # ? false


theta_Max = [
    [3.14159, 0.684719, 2.03444, 2.18628], # M=4
    [3.14155, 2.57765, 2.28452, 0.684719, 1.95839, 2.03444], # M=6
    [3.14159, 0.490883, 0.729713, 1.91063, 0.955239, 0.857072, 2.57765, 1.1832, 2.03444] # M=8
]
phi_Max = [
    [1.5708, 4.71239, 1.5708, 0.886077], # M=4
    [1.5708, 1.5708, 1.5708, 5.27633, 1.5708, 2.13474], # M=6
    [1.5708, 4.71239, 4.71199, 1.5708, 4.71282, 4.22151, 1.07991, 4.22151, 0.841069] # M=8
]

theta_Arc = [
    [0., 0.589437, 1.57448, 2.02376], # M=4
    [0., 0.384521, 0.927303, 0.135409, 1.57322, 1.55268], # M=6
    [3.14159, 2.85885, 0.674741, 0.0587125, 2.03432, 2.27237, 1.56898, 1.57232, 1.78567] # M=8
]

phi_Arc = [
    [1.57079, 4.71239, 1.5708, 6.22436], # M=4
    [1.5708, 4.71239, 4.71239, 6.13099, 1.5708, 3.35395], # M=6
    [1.5708, 1.5708, 1.5708, 2.84783, 1.5708, 5.20013, 4.71239, 1.90366, 0.0886607] # M=8    
]

theta_Det = [
    [1.78445, 1.28569, 0.00017959, 0.701871], # M=4
    [2.41051, 1.83873, 2.68454, 0.836022, 0.069125, 0.411515], # M=6
    [1.96918, 0.934956, 1.90402, 0.284836, 0.708879, 0.18032, 0.363484, 1.7756, 2.4393] # M=8
]

phi_Det = [
    [0.790252, 1.50861, 1.36399, 2.03546], # M=4
    [1.68016, 2.1397, 5.38404, 3.71373, 5.01612, 4.59139], # M=6
    [5.45757, 3.62961, 4.55925, 1.24255, 4.07992, 5.22654, 2.21233, 4.65093, 0.740643] # M=8
]

if angle_name == "Max"
    theta = theta_Max[M ÷ 2 - 1]
    phi = phi_Max[M ÷ 2 - 1]
elseif angle_name == "Arc"
    theta = theta_Arc[M ÷ 2 - 1]
    phi = phi_Arc[M ÷ 2 - 1]
elseif angle_name == "Det"
    theta = theta_Det[M ÷ 2 - 1]
    phi = phi_Det[M ÷ 2 - 1]
else
    error("Invalid angle_name: $angle_name. Must be one of \"Max\", \"Arc\", or \"Det\".")
end

equations = GramianMomentEquations1D3V(M, Kn, "ExtGram", slab_geometry=slab_geometry, theta=theta, phi=phi)

f_left = Maxwellian1D3D(7.0, (0.0, 0.0, 0.0), 1.0)
f_right = Maxwellian1D3D(1.0, (0.0, 0.0, 0.0), 1.0)

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
    variable = (u, eqns)->u[1]*u[5] # * Placeholder for the moment
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

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=cfl)

name = "1D3V/1D3V_angles/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_anglepair$(angle_name)"
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
