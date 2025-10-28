# ToDo: Implement testing for odd case

# miscellaneous testing stuff
test_closure(u, equations::GramianMomentEquations1D) = test_closure(u, Val(length(u)-1), Val(equations.closure))

# todo: update to new closure
# M3, extended
# function test_closure(u, M::Val{3}, closure::Val{:ExtGramOdd}) # Todo: <:Vector, <:SVector?? for u
#     @assert length(u) == 4
#     M = 3
#     n = (M+1)/2
#     χ = (n+1)/(2n)
#     A = u[4]*u[2]/u[1]
#     B = (u[4] - u[3]*u[2]/u[1])^2
#     C = u[3] - u[2]*u[2]/u[1]
#     return A + χ * B / C
# end

# M3, standard
function test_closure(u, M::Val{3}, closure::Val{:GramOdd})
    @assert length(u) == 4
    A = [u[3], u[4]]
    invG = inv([
    u[1] u[2];
    u[2] u[3]
    ])
    return A' * invG * A
    return A + χ * B / C
end

# M4, extended
function test_closure(u, M::Val{4}, closure::Val{:ExtGramEven})
    @assert length(u) == 5
    G = inv([u[1] u[2]; u[2] u[3]])
    a = [u[4], u[5]]'*G*[u[3], u[4]]
    b = (u[4]-u[3]*1/u[1]*u[2])/(u[3] - u[2]*1/u[1]*u[2])
    c = u[5] - [u[3], u[4]]'*G*[u[3], u[4]]
    return a + 1.5*b*c
end

# M4, standard
function test_closure(u, M::Val{4}, closure::Val{:GramEven})
    @assert length(u) == 5
    G = inv([u[1] u[2]; u[2] u[3]])
    a = [u[4], u[5]]'*G*[u[3], u[4]]
    return a
end

# todo: update to new closure
# M5, extended
# function test_closure(u, M::Val{5}, closure::Val{:ExtGramOdd})
#     @assert length(u) == 6
#     A_ = u[5]*(u[3]^2-u[2]*u[4])
#     B_ = u[6]*(-u[2]*u[3]+u[1]*u[4])
#     C_ = 1/(u[1]*u[3]-u[2]^2)
#     A = C_ * (A_ + B_)
#     D = u[5] - C_ * (u[3]*(u[3]^2 - u[2]*u[4]) + u[4] * (-u[2]*u[3] + u[1]*u[4]))
#     B = u[6] - C_ * (u[4]*(u[3]^2-u[2]*u[4]) + u[5]*(-u[2]*u[3]+u[1]*u[4]))
#     M = 5
#     n = (M+1)/2
#     χ = (n+1)/(2n)
#     return A + χ * B^2 / D
# end

function test_closure(u, M::Val{5}, closure::Val{:GramOdd})
    @assert length(u) == 6
    A = [u[4], u[5], u[6]]
    invG = inv([
    u[1] u[2] u[3];
    u[2] u[3] u[4];
    u[3] u[4] u[5]
    ])
    return A' * invG * A
end

function test_closure(u, M::Val{6}, closure::Val{:ExtGramEven})
    @assert length(u) == 7
    b = inv([
        u[1] u[2] u[3];
        u[2] u[3] u[4];
        u[3] u[4] u[5]
    ]) * [u[4]; u[5]; u[6]]
    c = inv([
        u[1] u[2];
        u[2] u[3]
    ]) * [u[3]; u[4]]
    A1 = [u[5]; u[6]; u[7]]'

    n = 3
    χ = (n+1)/n

    return A1 * b + χ * (
        u[6] - [u[4]; u[5]]' * c
    ) / (
        u[5] - [u[3]; u[4]]' * c
    ) * (
        u[7] - [u[4]; u[5]; u[6]]' * b
    )
end

function test_closure(u, M::Val{6}, closure::Val{:GramEven})
    @assert length(u) == 7
    b = inv([
        u[1] u[2] u[3];
        u[2] u[3] u[4];
        u[3] u[4] u[5]
    ]) * [u[4]; u[5]; u[6]]
    A1 = [u[5]; u[6]; u[7]]'

    return A1 * b
end


function check_realizability(u; verbose=true)
    M = length(u)-1
    if iseven(M)
        n = Int(M/2)
    else
        # ToDo: Can this be correct to manipulate n? Otherwise G_n is not defined
        n = Int((M+1)/2) - 1
    end
    # if verbose println("M: $M, n: $n, length(u): $(length(u))") end
    spd = all(eigen(gramian(u, n)).values .>= 0.0) #  all eigenvalues of gramian(u, n) should be positive (and real) for the matrix to be (semi-)positive definite
    if verbose println("Gramian positive semi-definite: ", spd) end

    if length(u) == 5
    cond1 = u[1] >= 0
    cond2 = u[1]*u[3] - u[2]^2 >= 0
    cond3 = u[1]*u[3]*u[5] - u[1]*u[4]^2 - u[2]^2*u[5] + 2*u[2]*u[3]*u[4] - u[3]^3 >= 0
    if verbose
        println("u[0] >= 0: ", cond1)
        println("u[0]*u[2] - u[1]^2 >= 0: ", cond2)
        println("u[0]*u[2]*u[4] - u[0]*u[3]^2 - u[1]^2*u[4] + 2*u[1]*u[2]*u[3] - u[2]^3 >= 0: ", cond3)
    end
    @assert (cond1 && cond2 && cond3) == spd # the three conditions should be equivalent to positive semi-definiteness
    end

    return spd
end




function foo(n=1_000_000, M=4)
    total = 0; valid = 0
    for i=1:n
        u = rand(M+1)
        ok = check_realizability(u, verbose=false)
        if ok
            total += 1
            J = flux_jacobian(u)
            if all(isreal.(eigen(J).values)) valid += 1 end
        end
    end
    return total, valid
end
