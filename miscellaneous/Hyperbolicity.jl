using Revise
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using ExtGram, Trixi
using LinearAlgebra

closure = "ExtGram"
ρ, v, θ = 1.0, 0.0, 1.0

Kn = 1.0 # Doesn't matter for the closure, necessary for defining the equations

# for M in 3:2:25
for M in 2:2:24
    equations = GramianMomentEquations1D(M, Kn, closure)
    for ρ in LinRange(0.1, 10.0, 11), v in LinRange(-5.0, 5.0, 11), θ in LinRange(0.1, 10.0, 11)
        f1 = Maxwellian(ρ, v, θ)
        f2 = Maxwellian(1.0, 0.0, 1.0)
        momentList = convective_moments(f1, Val(M+1)) .+ convective_moments(f2, Val(M+1))

        # Check realizability
        gramian = ExtGram.gramian(momentList, equations.n-1)
        # @assert isposdef(gramian) "Gramian is not positive definite for M = $M with moments $momentList"
        spd = all(eigen(ExtGram.gramian(momentList, equations.n)).values .>= 0.0)
        println("eigen(ExtGram.gramian(momentList, equations.n)).values = ", eigen(ExtGram.gramian(momentList, equations.n)).values)
        @assert spd "Gramian is not positive semi-definite for M = $M with moments $momentList"

        DF = ExtGram.flux_jacobian(momentList, equations)
        eigenvalues = eigvals(DF)

        number_distinct_eigenvalues = length(unique(eigenvalues))
        @assert number_distinct_eigenvalues == M + 1 "Expected $M + 1 distinct eigenvalues, got $number_distinct_eigenvalues"
        # println("Number of distinct eigenvalues for M = $M: $number_distinct_eigenvalues")
        real_eigenvalues = filter(x -> isreal(x), eigenvalues)
        # @assert length(real_eigenvalues) == M + 1 "Expected all eigenvalues to be real for M = $M, but got $(length(real_eigenvalues)) real eigenvalues with eigenvalues: $eigenvalues"
    # println("Number of real eigenvalues for M = $M: $(length(real_eigenvalues))")
        imag_eigenvalues = filter(x -> !isreal(x), eigenvalues)
        @assert length(imag_eigenvalues) == 0 "Expected no complex eigenvalues for M = $M, but got $(length(imag_eigenvalues)) with ρ=$(ρ) v=$(v), θ=$(θ) and eigenvalues: $eigenvalues"
    # println("Number of complex eigenvalues for M = $M: $(length(imag_eigenvalues))")
    end
end