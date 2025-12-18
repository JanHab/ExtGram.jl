"""
    Numerical utilities for 1D velocity - 3D moment transformations.

    The (extended) Gramian closure is constructed for the 1D1D system.
    However, we wish to extend this to a slab geometry in 3D moment space.

    Therefore, the moments are rotated.
    For M=4 we get the system:
    ```
    ∂_t U_{000} + ∂_x U_{100} &= S_{000}
    ∂_t U_{100} + ∂_x U_{200} &= S_{100}
    ∂_t U_{200} + ∂_x U_{300} &= S_{200}
    ∂_t U_{020} + ∂_x U_{120} &= S_{020}
    ∂_t U_{300} + ∂_x U_{400} &= S_{300}
    ∂_t U_{120} + ∂_x U_{220} &= S_{120}
    ∂_t U_{400} + ∂_x U_{500} &= S_{400}
    ∂_t U_{220} + ∂_x U_{320} &= S_{220}
    ∂_t U_{040} + ∂_x U_{140} &= S_{040}
    ∂_t U_{022} + ∂_x U_{122} &= S_{022}

    The moments U_{500}, U_{320}, U_{140}, U_{122} are not part of the slab geometry and calculated with the transformation.
"""

function idx(n::Integer)
    """
    Generates all sorted tuples of length 3 with elements in {1,2,3}
    """
    tuples = collect(Iterators.product(fill(1:3, n)...))
    sorted_vectors = [sort(collect(t)) for t in tuples]
    unique(sorted_vectors)
end

function index(n::Integer)
    """
    Gives the indices of the higher dimensional moments
    u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
    """
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

# Checks if the moment is non-zero in slab geometry
@inline slab_condition(y, z) = iseven(y) && iseven(z) && (y >= z)

function index_1d(n::Integer)
    """
    Positions of moments in slab geometry which are ≠ 0
    """
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
    """
    Generates the full list of multi-indices up to nmax
    """
    vcat([index(n) for n in 0:nmax]...)
end

function get_rotation_matrix(theta::Real, phi::Real)
    """
    Constructs the rotation matrix for given angles theta and phi.
    """
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
    """
    Constructs the tensor transformation matrix for rank-n tensors
    Combines rotation matrices for each index of the tensor.
    """
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
    """
    Helper function to construct a block diagonal matrix from a list of matrices.
    """
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
    """
    Constructs the full rotation matrix for moments up to given degree.
    """
    blocks = Matrix{Float64}[]
    push!(blocks, reshape([1.0], 1, 1)) # Degree 0
    
    for n in 1:moment_degree
        push!(blocks, tensor_transformation(n, theta, phi))
    end
    
    return construct_block_diagonal(blocks)
end

function nidx(moment_degree::Integer)
    """
    Index of moments necessary for the closure (u, u_x, u_xx, u_xxx, …)
    """ 
    Int[1 + div(k * (k + 1) * (k + 2), 6) for k in 0:moment_degree]
end

function shell_slab_positions(moment_degree::Integer)
    start = nidx(moment_degree)[end]
    [start + pos - 1 for pos in index_1d(moment_degree)]
end

function compile_fp(degree::Int)
    """
    Compiles the function to get the transformed higher order moments for given degree.
    """
    target_index = nidx(degree)[end]
    shell_positions = shell_slab_positions(degree)
    function fp(theta::Real, phi::Real)
        R = rot(degree, theta, phi)
        R[target_index, shell_positions]
    end
end

function get_valid_indices(M::Integer)
    """
    Get the valid indices for slab geometry up to moment degree M (which are ≠ 0)
    """
    slab_indices = Int[]
    index_start = 0
    for i in 0:M
        append!(slab_indices, index_1d(i) .+ index_start)
        index_start += length(index(i))
    end
    slab_indices
end
