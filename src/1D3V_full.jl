
function U_t_index(nmax::Integer)
    return mainmomindex(nmax)
end

function U_x_index(nmax::Integer)
    U_t = mainmomindex(nmax)
    return [U_t[i] .+ [1, 0, 0] for i in 1:length(U_t)]
end