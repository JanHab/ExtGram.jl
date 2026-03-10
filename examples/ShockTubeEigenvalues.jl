if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Revise, ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables, LinearAlgebra

# Parameter
M = 5    # number of moments
closure = "ExtGram" # flag for closure: "Gram", "ExtGram", "Grad"
Kn = 1.0  # Knudsen number
T_end = 0.3 
source = relaxation_source

# Domain and discretization parameters
x_lower = -2.0; x_upper = 2.0
domain = (x_lower, x_upper)
polydeg = 1  # polynomial degree
base_tree_level = 8 # initial mesh refinement level
surface_flux = flux_lax_friedrichs
volume_flux = flux_central

ρ_L = 7.0
v_L = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R = 0.0
θ_R = 1.0


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

# Store primitive variables in CSV file
CSV.write(
    "out/ShockTube/ρ_v_p.csv",
    Tables.columntable((x=x, rho=ρ, v=v, p=p))
)



# Eigenvalues
u_final = sol.u[end]
L = length(u_final)
Nloc = L ÷ (M+1)
Fmat = reshape(u_final, M+1, Nloc)

# x_index_center = Nloc ÷ 2
x_vector = [];
λ_vector = [];
u_vector = [];
for x_index in 1:Nloc
    # println("At x = $(x[x_index]):")
    u = SVector{M+1}(Fmat[:, x_index]...)
    jacobian = flux_jacobian(u, equations)
    λ = real.(eigen(jacobian).values)
    
    push!(x_vector, x[x_index])
    push!(λ_vector, λ...)
    push!(u_vector, u)
end

CSV.write(
    "out/ShockTube/eigenvalues.csv",
    Tables.columntable((x_vector=x_vector, eigenvalue=λ_vector, moments=u_vector))
)
