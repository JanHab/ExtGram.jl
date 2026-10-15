
using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi, OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays

closure = "ExtGram" # Doesn't matter for A
Kn = 1.0 # Doesn't matter for A
slab_geometry = false

M = 8

_U_t_index = ExtGram.U_t_index(M, slab_geometry)
_U_x_index = ExtGram.U_x_index(M, slab_geometry)
N_equations = length(_U_t_index)

position_in_u = Dict(index => i for (i, index) in enumerate(_U_t_index))
flux_index = Vector{Int}(undef, N_equations)
_closure_index = Vector{Vector{Int}}()
for i in 1:N_equations
    j = get(position_in_u, _U_x_index[i], 0)
    if j != 0
        flux_index[i] = j
    else
        push!(_closure_index, _U_x_index[i])
        flux_index[i] = -length(_closure_index)
    end
end
N_closure = length(_closure_index)

shell = slab_geometry ? ExtGram.index(M + 1)[ExtGram.index_1d(M + 1)] : ExtGram.index(M + 1) # ?index(M + 1)
N_shell = length(shell)
@assert _closure_index == shell[1:N_closure] "closure moments are not the leading block of index(M+1); the prefix assumption in closure_moments() is violated."

theta, phi = fibonacci_hemisphere_angles(N_shell)

N_full = length(ExtGram.mainmomindex(M))
S_state = slab_geometry ? ExtGram.slab_scatter_matrix(ExtGram.mainmomindex(M), _U_t_index) : Matrix{Float64}(I, N_full, N_full)
S_shell = slab_geometry ? ExtGram.slab_scatter_matrix(ExtGram.index(M + 1), shell) : Matrix{Float64}(I, N_shell, N_shell)

target_rows = ExtGram.nidx(M)
_P = [Matrix{Float64}(ExtGram.rot(M, theta[i], phi[i])[target_rows, :] * S_state) for i in 1:N_shell]

# Linear map from the order-(M+1) shell onto the directional closure values
T_shell = Matrix{Float64}(undef, N_shell, N_shell)
for i in 1:N_shell
    T_shell[i, :] = S_shell' * ExtGram.tensor_transformation(M + 1, theta[i], phi[i])[1, :]
end

# @assert (all(sum(abs.(T_shell),dims=2) .<= 2abs.(diag(T_shell))))
heatmap(1:N_shell, 1:N_shell, T_shell)