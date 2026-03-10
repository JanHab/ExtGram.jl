"""
    Example file for a 1D-1D shock (tube or structure) test case with the DVM method.
"""

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
N = parse(Int, ARGS[1])
c_l = parse(Float64, ARGS[2])
c_u = parse(Float64, ARGS[3])
Kn = parse(Float64, ARGS[4])
source_string = ARGS[5]
source = source_string == "relaxation_source" ? relaxation_source : zero_source 
T_end = parse(Float64, ARGS[6])
base_tree_level = parse(Int, ARGS[7]) # e.g. 10
polydeg = parse(Int, ARGS[8]) # e.g. 1
ρ_L = parse(Float64, ARGS[9]) # 7.0
v_L = parse(Float64, ARGS[10]) # 0.0
θ_L = parse(Float64, ARGS[11]) # 1.0
ρ_R = parse(Float64, ARGS[12]) # 1.0
v_R = parse(Float64, ARGS[13]) # 0.0
θ_R = parse(Float64, ARGS[14]) # 1.0
x_left = parse(Float64, ARGS[15]) # -5.0
x_right = parse(Float64, ARGS[16]) # 5.0
domain = (x_left, x_right)

# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupDVM1DShockTube(
    N, Kn,
    c_l, c_u, 
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R);
    base_tree_level = base_tree_level,
    domain = domain,
    polydeg = polydeg,
)

semi = SemidiscretizationHyperbolic(
    mesh, equations, 
    initial_condition, solver, 
    boundary_conditions=boundary_conditions, 
    source_terms=source
)

#= set up ODE =#
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)

# Some parameters for callbacks
cfl = 0.99          # Maximum cfl number
plot_interval = 20  # plot every 20 steps
time_interval = 20
name="ShockTube/DVM/DVM_solution_N$(N)_c_l$(c_l)_c_u$(c_u)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_x_left$(x_left)_x_right$(x_right)" # name of output files

alive_callback = AliveCallback(analysis_interval=100)
summary_callback = SummaryCallback()
stepsize_callback = StepsizeCallback(cfl=cfl)

save_solution = SaveTriangulationCallback(
    time_interval=tspan[2]/time_interval,
    save_initial_solution=true,
    file_format="tsv",
    append_solution=true,
    solution_variables = cons2cons,
    clear_out_dir=false,
    name=name,
    info="basis = $(Base.typename(typeof(basis)).wrapper)"
)

callbacks = CallbackSet(
    alive_callback,
    stepsize_callback,
    save_solution,
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

x, ρ, v, θ, p, q = ρ_v_θ_p_DVM(semi, sol, equations)

# Store primitive variables in CSV file
CSV.write(
    "out/ShockTube/DVM/DVM_solution_N$(N)_c_l$(c_l)_c_u$(c_u)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_x_left$(x_left)_x_right$(x_right)_rho_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p, q=q))
)