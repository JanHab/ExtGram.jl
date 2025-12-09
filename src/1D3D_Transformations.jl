# Numerical utilities for 1D velocity / 3D moment transformations.

function idx(n::Integer)
    tuples = collect(Iterators.product(fill(1:3, n)...))
    sorted_vectors = [sort(collect(t)) for t in tuples]
    unique(sorted_vectors)
end

function index(n::Integer)
    tuples = collect(Iterators.product(0:n, 0:n, 0:n))
    vectors = [collect(t) for t in tuples]
    selected = filter(v -> sum(v) == n, vectors)
    sort!(selected)
    reverse(selected)
end

@inline function tuple_to_powers(tup)
    c1 = count(==(1), tup)
    c2 = count(==(2), tup)
    c3 = count(==(3), tup)
    [c1, c2, c3]
end

@inline slab_condition(y, z) = iseven(y) && iseven(z) && (y >= z)

function index_1d(n::Integer)
    indices = index(n)
    pos = Int[]
    for (idx_pos, (_, y, z)) in enumerate(indices)
        if slab_condition(y, z)
            push!(pos, idx_pos)
        end
    end
    pos
end

function mainmomindex(nmax::Integer)
    vcat([index(n) for n in 0:nmax]...)
end

function get_rotation_matrix(theta::Real, phi::Real)
    θ = Float64(theta)
    φ = Float64(phi)
    ct, st = cos(θ), sin(θ)
    cp, sp = cos(φ), sin(φ)
    [
        -ct * sp   -st * sp    cp;
         st * sp   -ct * sp    0.0;
        -ct * cp   -cp * st   -sp
    ]
end

function tensor_transformation(n::Int, theta::Real, phi::Real)
    row_basis = idx(n) 
    col_basis = index(n) 
    col_lookup = Dict(vec => i for (i, vec) in enumerate(col_basis))
    
    dim_rows = length(row_basis)
    dim_cols = length(col_basis)
    
    # Initialize with Float64 zeros instead of generic/symbolic
    TT = zeros(Float64, dim_rows, dim_cols)
    R = get_rotation_matrix(theta, phi)
    
    # Pre-allocate tuples generator to avoid recreating it inside loops if possible
    all_tuples = collect(Iterators.product(fill(1:3, n)...))

    for (r_idx, iList) in enumerate(row_basis)
        for jList_tuple in all_tuples
            # Compute Coefficient
            coef = 1.0
            for k in 1:n
                # iList is the row index (target), jList is the column index (source)
                coef *= R[iList[k], jList_tuple[k]]
            end
            
            # Identify Column
            # tuple_to_powers inline logic for speed
            c1 = count(==(1), jList_tuple)
            c2 = count(==(2), jList_tuple)
            c3 = count(==(3), jList_tuple)
            powers = [c1, c2, c3]
            
            if haskey(col_lookup, powers)
                c_idx = col_lookup[powers]
                TT[r_idx, c_idx] += coef
            end
        end
    end
    return TT
end

function construct_block_diagonal(matrices::Vector{Matrix{Float64}})
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)
    
    full_matrix = zeros(Float64, total_rows, total_cols)
    
    current_row = 1
    current_col = 1
    
    for m in matrices
        r, c = size(m)
        full_matrix[current_row:(current_row + r - 1), 
                    current_col:(current_col + c - 1)] = m
        current_row += r
        current_col += c
    end
    return full_matrix
end

function rot(moment_degree::Int, theta::Real, phi::Real)
    blocks = Matrix{Float64}[]
    push!(blocks, reshape([1.0], 1, 1)) # Degree 0
    
    for n in 1:moment_degree
        push!(blocks, tensor_transformation(n, theta, phi))
    end
    
    return construct_block_diagonal(blocks)
end

function nidx(moment_degree::Integer)
    Int[1 + div(k * (k + 1) * (k + 2), 6) for k in 0:moment_degree]
end

function shell_slab_positions(moment_degree::Integer)
    start = nidx(moment_degree)[end]
    [start + pos - 1 for pos in index_1d(moment_degree)]
end

function compile_fp(degree::Int)
    target_index = nidx(degree)[end]
    shell_positions = shell_slab_positions(degree)
    function fp(theta::Real, phi::Real)
        R = rot(degree, theta, phi)
        R[target_index, shell_positions]
    end
end

function get_valid_indices(M::Integer)
    slab_indices = Int[]
    index_start = 0
    for i in 0:M
        append!(slab_indices, index_1d(i) .+ index_start)
        index_start += length(index(i))
    end
    slab_indices
end
