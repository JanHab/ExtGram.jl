"""
    Example file for a 1D-1D shock (tube or structure) test case with the moment method.
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
source = source_string == "relaxation_source" ? relaxation_source : zero_source # default to zero_source if not relaxation_source
T_end = parse(Float64, ARGS[5])
base_tree_level = parse(Int, ARGS[6]) # e.g. 10
polydeg = parse(Int, ARGS[7]) # e.g. 1
ρ_L = parse(Float64, ARGS[8]) # 7.0
v_L = parse(Float64, ARGS[9]) # 0.0
θ_L = parse(Float64, ARGS[10]) # 1.0
ρ_R = parse(Float64, ARGS[11]) # 1.0
v_R = parse(Float64, ARGS[12]) # 0.0
θ_R = parse(Float64, ARGS[13]) # 1.0
x_lower = parse(Float64, ARGS[14]) # -2.0
x_upper = parse(Float64, ARGS[15]) # 2.0

# Fixed settings
domain = (x_lower, x_upper)



# Setting up everything
basis, mesh, equations, initial_condition, solver, boundary_conditions = setupGramianMomentEquations1DShockTube(
    M, Kn, closure,
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

callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.45,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="ShockTube/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)", # name of output files
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

# Eigenvalues
if polydeg == 1 # only implemented for polydeg=1 (linear basis functions)
    u_final = sol.u[end]
    L = length(u_final)
    Nloc = L ÷ (M+1)
    Fmat = reshape(u_final, M+1, Nloc)

    x_vector = [];
    λ_vector = [];
    u_vector = [];
    for x_index in 1:Nloc
        u = SVector{M+1}(Fmat[:, x_index]...)
        jacobian = flux_jacobian(u, equations)
        λ = real.(eigen(jacobian).values)
        
        push!(x_vector, x[x_index])
        push!(λ_vector, λ)
        push!(u_vector, u)
    end

    CSV.write(
        "out/ShockTube/Moments/gram_solution_M$(M)_closure$(closure)_Kn$(Kn)_source$(source_string)_T_end$(T_end)_rho_L$(ρ_L)_rho_R$(ρ_R)_v_L$(v_L)_v_R$(v_R)_theta_L$(θ_L)_theta_R$(θ_R)_base_tree_level$(base_tree_level)_polydeg$(polydeg)_eigenvalues.csv",
        Tables.columntable((x_vector=x_vector, eigenvalue=λ_vector, moments=u_vector))
    )
end