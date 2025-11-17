using Plots, Trapz

f0(v, ϕ, v0) = 1 / √(2π) * exp(-.5 * (sign(v) * √(v^2 - 2ϕ) - v0)^2)
f1(v, ϕ, v0, β) = 1 / √(2π) * exp(-.5 * (β  * (v^2 - 2ϕ) + v0)^2)
f(v, ϕ, v0, β) = if v^2 > 2ϕ; f0(v, ϕ, v0); else; f1(v, ϕ, v0, β); end;

# Parameter
ψ = 0.2
Δ = 500
x0 = 512
ϕ(x) = ψ * exp(((x-x0)/Δ)^2) #10.
v0 = 1.
β = -0.1

x = LinRange(-1000, 2000, 5_000)
v = LinRange(-20, 20, 1_000)

heatmap(x, v, [f(vj, ϕ(xi), v0, β) for xi in x, vj in v]', xlabel="x", ylabel="v", title="Electron hole distribution function f(x,v)", colorbar_title="f(x,v)")
# Z = [f(vj, ϕ(xi), v0, β) for xi in x, vj in v]'
# surface(x, v, Z,
#     xlabel="x", ylabel="v", zlabel="f(x,v)",
#     title="Electron hole distribution function f(x,v)",
#     colorbar=true)

ρ, v_mean, θ = zeros(length(x)), zeros(length(x)), zeros(length(x))
for x_ in x
    ρ_ = trapz(v, f.(v, ϕ(x_), v0, β))
    v_mean_ = trapz(v, v .* f.(v, ϕ(x_), v0, β)) / ρ_
    θ_ = trapz(v, (v .- v_mean_).^2 .* f.(v, ϕ(x_), v0, β)) / ρ_
    ρ[Int(round((x_ + 1000) / 3000 * (length(x)-1) + 1))] = ρ_
    v_mean[Int(round((x_ + 1000) / 3000 * (length(x)-1) + 1))] = v_mean_
    θ[Int(round((x_ + 1000) / 3000 * (length(x)-1) + 1))] = θ_
end

plot(x, ρ, xlabel="x", ylabel="Density ρ", title="Density profile of electron hole", legend=false)
plot(x, v_mean, xlabel="x", ylabel="Mean velocity v", title="Mean velocity profile of electron hole", legend=false)
plot(x, θ, xlabel="x", ylabel="Temperature θ", title="Temperature profile of electron hole", legend=false)

ϕ_ = 10.
v0 = 1.
β = -0.1

v = LinRange(-10, 10, 500)
plot(v, f.(v, ϕ_, v0, β), xlabel="v", ylabel="f(v)", title="Electron hole distribution function", legend=false)

# Compute moments
ρ = trapz(v, f.(v, ϕ, v0, β))
v_mean = trapz(v, v .* f.(v, ϕ, v0, β)) / ρ
θ = trapz(v, (v .- v_mean).^2 .* f.(v, ϕ, v0, β)) / ρ

println("Density ρ = $ρ")
println("Mean velocity v = $v_mean")
println("Temperature θ = $θ")