"""
    Example file for a 1D-1D shock (tube or structure) test case,where the gauge parameter χ is specified via user input.
"""

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = parse(Int, ARGS[1])
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = parse(Float64, ARGS[3])
source_string = ARGS[4]
source = source_string == "relaxation_source" ? relaxation_source : zero_source 
T_end = parse(Float64, ARGS[5])
base_tree_level = parse(Int, ARGS[6]) # e.g. 10
polydeg = parse(Int, ARGS[7]) # e.g. 3
χ_value = ARGS[8] # chi value
χ = χ_value
if χ_value != "optimal"
    χ = parse(Float64, χ_value)
end
ρ_L = parse(Float64, ARGS[9]) # 7.0
v_L = parse(Float64, ARGS[10]) # 0.0
θ_L = parse(Float64, ARGS[11]) # 1.0
ρ_R = parse(Float64, ARGS[12]) # 1.0
v_R = parse(Float64, ARGS[13]) # 0.0
θ_R = parse(Float64, ARGS[14]) # 1.0

# Fixed settings
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)


# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DShockTube(
    M, Kn, closure,
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R);
    base_tree_level = base_tree_level,
    polydeg = polydeg,
    domain = domain,
    χ_set = χ # set χ value explicitly
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

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="Shock/GaugeInvestigation/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_chi$(χ_value)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name of output files
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
