using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Access arguments by index
M = parse(Int, ARGS[1])
closure = ARGS[2] # String # "Gram", "ExtGram" or "Grad"
Kn = parse(Float64, ARGS[3])
source_string = ARGS[4]
source = relaxation_source
T_end = parse(Float64, ARGS[5])
base_tree_level = parse(Int, ARGS[6]) # e.g. 8
polydeg = parse(Int, ARGS[7]) # e.g. 3
ρ_L = parse(Float64, ARGS[8]) # 7.0
v_L = parse(Float64, ARGS[9]) # 0.0
θ_L = parse(Float64, ARGS[10]) # 1.0
ρ_R = parse(Float64, ARGS[11]) # 1.0
v_R = parse(Float64, ARGS[12]) # 0.0
θ_R = parse(Float64, ARGS[13]) # 1.0


# Fixed settings
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)



# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DRiemann(
    M, Kn, closure,
    Maxwellian(ρ_L, v_L, θ_L), # Density, velocity, temperature
    Maxwellian(ρ_R, v_R, θ_R);
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

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution"
    # name="Riemann1D/Moments/relaxation/gram_solution_M$(M)_Kn$(Kn)_source$(ARGS[4])_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name of output files
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
savefig(p1, "out/Riemann1D/Moments/relaxation/gram_solution_M$(M)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_rho_v_p.pdf")

# Store primitive variables in CSV file
CSV.write(
    "out/Riemann1D/Moments/relaxation/gram_solution_M$(M)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_rho_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)

# Plot maximum eigenvalue (wave-speed) of flux Jacobian over time
n_plots = 5
p2 = plot_λ_max(semi, sol, M, n_plots, x_lower, x_upper)
display(p2)
savefig(p2, "out/Riemann1D/Moments/relaxation/gram_solution_M$(M)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_lambda_max.pdf")