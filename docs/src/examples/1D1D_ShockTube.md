# Shock Tube Example

This tutorial mirrors the setup from examples/Shock_1D1D.jl with parameters chosen for the shock tube test case. 
It evolves a one-dimensional shock tube and stores the resulting density/velocity/pressure profile.
Note: The grid is usually finer, and is coarsened for this example.

``` julia
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, CSV, Tables
```

### Setup

Start with the number of moments, closure, and the Knudsen number

``` julia
M = 4
closure = "ExtGram" # possible: "Gram", "ExtGram" or "Grad"
Kn = 1.0
source = relaxation_source # if "zero_source" Kn doesn't matter, as no source term
```

continue with the numerical setting
``` julia
T_end = 0.3
base_tree_level = 8 # 10 used in thesis, but let's speed it up a bit
polydeg = 1
x_lower = -2.0
x_upper = 2.0

# Fixed settings
domain = (x_lower, x_upper)
```

and the initial conditions
```julia
ρ_L = 7.0
v_L = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R = 0.0
θ_R = 1.0
```

### Set up the numerical solver

``` julia
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
```

with the ode problem

``` julia
tspan = (0.0, T_end)
ode = semidiscretize(semi, tspan)
```

and the callbacks

``` julia
callbacks, summary_callback = callbacksGramianMomentEquations(
    semi, tspan, basis; 
    cfl = 0.9,          # Maximum cfl number
    plot_interval = 20,  # plot every 20 steps
    name="gram_solution",
)
```

### Solve the problem

``` julia
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
```

### Post Processing

Visualize the density, velocity, and pressure

``` julia
x, ρ, v, p, p1 = plot_ρ_v_p(sol, M, x_lower, x_upper)
plot!(
    p1, xlims = (-1, 1)
)
display(p1)
```

![ShockTube](assets/1D1D_ShockTube.svg)

or visualize all the (conservative) moments separately

``` julia
u_final = sol.u[end]
L = length(u_final)
Nloc = L ÷ (M + 1)
Fmat = reshape(u_final, M+1, Nloc)

u = []
for i in 1:(M+1)
    push!(u, Fmat[i, :])
end

# Plot first moments
pl = Vector{Any}(undef, M+1)
x = range(x_lower, x_upper, length=Nloc)

for i in 1:(M+1)
    pl[i] = scatter()
    plot!(
        pl[i],
        xlabel="x",
        ylabel="u^{$(i-1)}",
    )
    plot!(
        pl[i],
        x, u[i], 
        label="U$(i-1)",
    )
end

display(pl)
l = @layout  [grid(3,3)]
plot(pl..., layout=l, size=(1_200, 600))
```

![ShockTube](assets/1D1D_ShockTube_Moments.svg)