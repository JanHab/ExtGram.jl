using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi

M = 4
closure = "ExtGram"
Kn = 1.0
source = relaxation_source

# IC
ρ_L = 7.0
v_L1 = 0.0
v_L2 = 0.0
v_L3 = 0.0
θ_L = 1.0
ρ_R = 1.0
v_R1 = 0.0
v_R2 = 0.0
v_R3 = 0.0
θ_R = 1.0

equations = GramianMomentEquations1D3D(M, Kn, closure)

initial_condition = InitialConditionsShockTube1D3D(
    Maxwellian1D3D(ρ_L, (v_L1, v_L2, v_L3), θ_L), # Density, velocity, temperature
    Maxwellian1D3D(ρ_R, (v_R1, v_R2, v_R3), θ_R), # Shock in density, but not velocity, temperature initially
    M,
    equations
)

left_full = convective_moments_1D3D(M, ρ_L, (v_L1, v_L2, v_L3), θ_L)
slab_indices = HyQMOM.get_valid_indices(M)
u_cons = SVector{length(slab_indices)}(left_full[slab_indices])

w_prim = cons2prim(u_cons, equations)

# closure_transformation_prim = HyQMOM.closure_transform(w_prim, equations)
closure_transformation_cons = HyQMOM.closure_transform(u_cons, equations)
extended_moments_cons = SVector{length(slab_indices)+length(closure_transformation_cons)}(u_cons..., closure_transformation_cons...)

extended_moments_prim = HyQMOM.moment_cons2prim(extended_moments_cons, equations, MOMENT_INDICES_CLOSURE)

# extended_moments_prim = SVector{14}(w_prim..., closure_transformation_prim...)

MOMENT_INDICES_CLOSURE = (
    (0,0,0), # 1: rho
    (1,0,0), # 2: rho * v
    (2,0,0), # 3: P_xx
    (0,2,0), # 4: P_yy
    (3,0,0), # 5
    (1,2,0), # 6
    (4,0,0), # 7
    (2,2,0), # 8
    (0,4,0), # 9
    (0,2,2),  # 10: P_yyzz (mixed transverse)
    (5,0,0), # 11
    (3,2,0), # 12
    (1,4,0), # 13
    (1,2,2)  # 14
)

# T = eltype(extended_moments_prim)
# v = w_prim[2]
# closure_transformation_cons = SVector{14, T}(ntuple(idx -> begin
#     i, j, k = MOMENT_INDICES_CLOSURE[idx]
#     val = zero(T)
#     for m in 0:i
#         # Look up P_{mjk} from the central moment vector
#         p_mjk = HyQMOM.get_moment_val(extended_moments_prim, m, j, k)
#         val += binomial(i, m) * v^(i-m) * p_mjk
#     end
#     return val
# end, 14))

# known_moments_cons = SVector{6}(ntuple(i->if i <= 2; u_cons[i+1]; else; u_cons[i+2]; end, 6))
# flux_moments_cons = SVector{10}(known_moments_cons..., closure_transformation_cons[11:14]...)

# u_cons
# flux_moments_cons