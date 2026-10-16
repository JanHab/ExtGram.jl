# Execute with "julia -t auto --project=. miscellaneous/DVM_ReferenceKinetics_jl.jl" to parallelize

# 1D-3V BGK shock tube with KitBase.jl (the core of Kinetic.jl),
using KitBase
using Plots
using DelimitedFiles

# ---------------- case parameters ----------------
const ρL   = 7.0        # 7.0 or 30.0
const ρR   = 1.0
const θ0   = 1.0        # temperature on both sides, θ = T
const Kn   = parse(Float64, ARGS[1])        # 0.2 or 1.0
const tend = 0.3
const x0, x1, nx = -2.0, 2.0, 1_024     # double nx once to check grid convergence
const umax, nu    = 10.0, 400          # double nu once to check velocity-grid convergence

# ---- setup: one distribution function, 3D velocity space ----
set = Setup(matter = "gas", case = "shocktube",
            space = "1d1f3v", flux = "kfvs", collision = "bgk", nSpecies = 1,
            interpOrder = 2, limiter = "vanleer", boundary = "fix",
            cfl = 0.5, maxTime = tend)

# K = 0: all three translational directions are resolved, nothing is "internal"; γ stays 5/3
gas = Gas(Kn = Kn, Ma = 0.0, Pr = 1.0, K = 0.0, γ = 5/3,
          ω = 1.0, αᵣ = 1.0, ωᵣ = 0.5, μᵣ = Kn / 2)

# ---- velocity space ----
# c1 needs the wide, fine grid; c2, c3 carry no drift and θ ≤ 1.5, so ±6 with ~32-48 points suffices
ps = PSpace1D(x0, x1, nx, 1)
vs = VSpace3D(-umax, umax, nu, -6.0, 6.0, 64, -6.0, 6.0, 64)

# ---- initial condition: prim = [ρ, v1, v2, v3, λ] ----
primL = [ρL, 0.0, 0.0, 0.0, 1 / (2θ0)]
primR = [ρR, 0.0, 0.0, 0.0, 1 / (2θ0)]
p = (wL = prim_conserve(primL, gas.γ), wR = prim_conserve(primR, gas.γ),
     primL = primL, primR = primR, γ = gas.γ, u = vs.u, v = vs.v, w = vs.w)

fw(x, p) = x < 0.0 ? p.wL   : p.wR          # (ρ, ρv1, ρv2, ρv3, ρE)
bc(x, p) = x < 0.0 ? p.primL : p.primR
ff(x, p) = maxwellian(p.u, p.v, p.w, conserve_prim(fw(x, p), p.γ))   # 3D array
ib = IB1F(fw, ff, bc, p)

ks = SolverSet(set, ps, vs, gas, ib)

# sanity check before trusting anything: must print your τ at ρ = 1, θ = 1
println("τ(ρ=1, θ=1) = ", vhs_collision_time(primR, ks.gas.μᵣ, ks.gas.ω))

# ---------------- time integration ----------------
function run!(ks)
    ctr, face = init_fvm(ks)            # older versions: init_fvm(ks, ks.ps)
    t, res = 0.0, zeros(3)
    while t < ks.set.maxTime
        dt = min(timestep(ks, ctr, t), ks.set.maxTime - t)   # land exactly on tend
        reconstruct!(ks, ctr)           # slopes for interpOrder = 2
        evolve!(ks, ctr, face, dt)
        update!(ks, ctr, face, dt, res)
        t += dt
    end
    return ctr
end
KitBase.vhs_collision_time(prim::AbstractVector{T}, μᵣ, ω) where {T<:Real} = Kn
println("τ(left) = ", vhs_collision_time(primL, ks.gas.μᵣ, ks.gas.ω), "   τ(right) = ", vhs_collision_time(primR, ks.gas.μᵣ, ks.gas.ω))
methods(KitBase.vhs_collision_time)
ctr = run!(ks)

# ---------------- output ----------------
u, v, w, wt = ks.vs.u, ks.vs.v, ks.vs.w, ks.vs.weights     # all 3D arrays
mom(a, b, c, f) = sum(@. u^a * v^b * w^c * f * wt)          # U_{abc}
out = zeros(ks.ps.nx, 9)
for i in 1:ks.ps.nx
    ρ, v1, v2, v3, λ = ctr[i].prim
    θ = 1 / (2λ)
    f = ctr[i].f
    # e.g. U500, U320, U140, U122 for the M = 4 closure target
    out[i, :] = [ks.ps.x[i], ρ, v1, θ, ρ * θ,
                 mom(5,0,0,f), mom(3,2,0,f), mom(1,4,0,f), mom(1,2,2,f)]
end
writedlm("bgk_3V_rhoL$(ρL)_Kn$(Kn).csv", out, ',')

plot(out[:,1], out[:,2], label="ρ", xlabel="x", ylabel="ρ", title="Shock Tube BGK 1D-3V Slab")
plot(out[:,1], out[:,3], label="v", xlabel="x", ylabel="v", title="Shock Tube BGK 1D-3V Slab")
plot(out[:,1], out[:,4], label="θ", xlabel="x", ylabel="θ", title="Shock Tube BGK 1D-3V Slab")
plot(out[:,1], out[:,5], label="ρθ", xlabel="x", ylabel="ρθ", title="Shock Tube BGK 1D-3V Slab")
