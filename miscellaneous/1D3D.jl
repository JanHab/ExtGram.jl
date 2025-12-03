# Disclaimer: This Julia code is a translation of a Mathematica script 
# The Translation was mainly done with Google Gemini and manually refined.
# Some constructs may not have direct equivalents and were adapted accordingly.

if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using Symbolics
using LinearAlgebra
using Distributions
using HyQMOM, Trixi

# --- 1. Index Generation Functions ---

# Equivalent to: idx[n_]
# Generates all sorted tuples of length 3 with elements in {1,2,3}
function idx(n)
    # 1. Generate tuples
    raw_tuples = vec(collect(Iterators.product(fill(1:3, n)...)))
    
    # Convert tuples to vectors and sort the elements inside each vector
    # e.g., (3, 1, 2) -> [1, 2, 3]
    sorted_vectors = [sort(collect(t)) for t in raw_tuples]
    
    # Remove duplicates
    return unique(sorted_vectors)
end
idx(4) # Test call

# Equivalent to: index[n_]
# Gives the indizes of the higher dimensional moments
# u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
function index(n)
    # Generate tuples 0 to n of length 3
    # Iterators.product creates the Cartesian product
    raw_tuples = collect(Iterators.product(0:n, 0:n, 0:n))
    vecs = [collect(t) for t in raw_tuples]
    selected = filter(v -> sum(v) == n, vecs) # only consider sum is n
    # Mathematica sorts lexicographically then Reverses. 
    # We sort standardly then reverse to match the script perfectly.
    sort!(selected)
    return reverse(selected)
end
index(4) # Test call

# ToDo: Is this necessary? What is it doing?
# Equivalent to: tupleToPowers[tup_, n_]
function tuple_to_powers(tup, n)
    c1 = count(==(1), tup)
    c2 = count(==(2), tup)
    c3 = count(==(3), tup)
    return [c1, c2, c3]
end

# Equivalent to: index1D[n_]
# Positions of moments in slab geometry which are ≠ 0
function index_1d(n)
    full_indices = index(n)

    # Only consider indices with y and z even and y >= z (no duplicates)
    pos = Int[]
    for (i, val) in enumerate(full_indices)
        _, y, z = val
        if iseven(y) && iseven(z) && (y >= z)
            push!(pos, i)
        end
    end
    return pos
end
index_1d(2) # Test call

# Equivalent to: mainmomindex[nmax_]
# Generates the full list of multi-indices up to nmax
function mainmomindex(nmax)
    return vcat([index(n) for n in 0:nmax]...)
end
mainmomindex(4) # Test call



# --- 3. Tensor Transformation Logic ---

# Memoization Cache: Key is (degree, theta_symbol, phi_symbol)
# ToDo: Do it completely numeric later
# ToDo: Work with cache?
# Constructs the rotation matrix for given theta, phi
function get_rotation_matrix(theta, phi)
    # If theta/phi are symbolic, these return symbolic expressions.
    # ToDo: Ensure compatibility with numeric inputs later.
    ct, st = cos(theta), sin(theta)
    cp, sp = cos(phi), sin(phi)
    
    return [
        -ct*sp  -st*sp   cp;
         st*sp  -ct*sp   0;
        -ct*cp  -cp*st  -sp
    ]
end
# Define symbolic variables
@variables θ φ 
get_rotation_matrix(θ, φ) # Test call

function tensor_transformation(n, theta, phi)
    row_basis = idx(n) # indices
    col_basis = index(n) # moments
    col_lookup = Dict(vec => i for (i, vec) in enumerate(col_basis)) # allocating indices to moments
    
    # CRITICAL CHANGE: Initialize with symbolic zeros (Num(0))
    # We use size matching row/col basis
    dim_rows = length(row_basis)
    dim_cols = length(col_basis)
    TT = zeros(Num, dim_rows, dim_cols)
    
    R = get_rotation_matrix(theta, phi)
    all_tuples = collect(Iterators.product(fill(1:3, n)...))

    # This is where the magic happens
    for (r_idx, iList) in enumerate(row_basis)
        for jList_tuple in all_tuples
            jList = collect(jList_tuple)
            
            # Calculate coefficient product
            # Initialize coef as generic 1 (integer) so it can become symbolic
            coef = 1
            for k in 1:n
                coef *= R[iList[k], jList[k]]
            end
            
            powers = tuple_to_powers(jList, n)
            
            if haskey(col_lookup, powers)
                c_idx = col_lookup[powers]
                TT[r_idx, c_idx] += coef
            end
        end
    end

    return TT
end
tensor_transformation(2, θ, φ) # Test call with symbolic angles
tensor_transformation(2, 0.5, 0.3) # Test call with numeric angles




# --- 4. Symbolic Block Diagonal Construction ---

# ToDo: Work with cache?

# Helper function to manually construct a Block Diagonal Matrix
# (Stitches a list of matrices [M1, M2...] into one large matrix)
function construct_block_diagonal(matrices::Vector{Matrix{Num}})
    # Calculate total dimensions
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)
    
    # Initialize large symbolic matrix with zeros
    # ToDo: Can assign Float64 later?
    full_matrix = zeros(Num, total_rows, total_cols) # 'Num' to hold symbolic variables
    
    current_row = 1
    current_col = 1
    
    for m in matrices
        rows, cols = size(m)
        # Paste the current matrix into the correct block position
        full_matrix[current_row:(current_row + rows - 1), 
                    current_col:(current_col + cols - 1)] = m
        
        # Move the cursor
        current_row += rows
        current_col += cols
    end
    
    return full_matrix
end

# Equivalent to: RotSymbolic[MomentDegree_]
# Generates the full block diagonal rotation matrix symbolically/numerically
function rot(moment_degree, theta, phi)

    # Collect the blocks for n = 0 to moment_degree
    blocks = Matrix{Num}[] # "Num" to hold symbolic entries, # ToDo: Can assign Float64 later?
    
    # First block matrix
    # IdentityMatrix[1] is just a 1x1 matrix containing 1
    push!(blocks, reshape([Num(1)], 1, 1)) # "Num" again
    for n in 1:moment_degree
        # Call the tensor transformation function we defined earlier
        mat = tensor_transformation(n, theta, phi)
        push!(blocks, mat)
    end
    
    # Stitch them together
    final_matrix = construct_block_diagonal(blocks)
    
    return final_matrix
end
rot(2, θ, φ) # Test call
rot(2, 0.5, 0.3) # Test call



# --- 2. VarsSymbolic Implementation ---
# Equivalent to: VarsSymbolic[MomentDegree]
# This generates the vector of variables V_i_j_k, which are ≠ 0 in slab geometry
function vars_symbolic(moment_degree)
    # Initialize an empty vector of Symbolic Numbers
    vars = Num[]
    
    # n=0 is always valid (V_0_0_0)
    push!(vars, first(@variables V_0_0_0))
    for n in 1:moment_degree
        # 1. Get the FULL list of indices for this degree
        full_indices = index(n) 
        
        for (i, j, k) in full_indices
            # 2. Check the validity condition:
            # EvenQ[y] && EvenQ[z] && (y >= z), Here: y->j, z->k
            if iseven(j) && iseven(k) && (j >= k)
                # Case A: Valid -> Generate Symbol
                sym_name = Symbol("V_$(i)_$(j)_$(k)")
                var = Symbolics.unwrap(first(@variables $sym_name))
                push!(vars, Num(var))
            else
                # Case B: Invalid -> Insert Zero
                push!(vars, Num(0))
            end
        end
    end
    
    return vars
end
vars_symbolic(2) # Test call




# Index of moments necessary for the closure (u, u_x, u_xx, u_xxx, …)
function nidx(moment_degree)
    # This generates the start index of every moment shell in the flattened list
    return Int[1 + div(k*(k+1)*(k+2), 6) for k in 0:moment_degree]
end
nidx(2) # Test call

# --- 2. Unknowns Helper (x) ---
# Which moments we are looking for in the closure
# ToDo: Pretty sure this is not necessary at all, but let's see
function get_unknowns_shell(n)
    vars = Num[]
    
    # 1. Iterate only indices for degree 'n'
    full_indices = index(n)
    
    for (i, j, k) in full_indices
        # 2. Check symmetry condition
        if iseven(j) && iseven(k) && (j >= k)
            # Create symbol
            sym_name = Symbol("V_$(i)_$(j)_$(k)")
            var = Symbolics.unwrap(first(@variables $sym_name))
            push!(vars, Num(var))
        end
    end
    
    return vars
end
get_unknowns_shell(2) # Test call



# This constructs a function to be evaluated on θ, φ for a given moment degree
# Returns the coefficients matrix
# * I am sure that this is not the most efficient way to do this
# --- 3. FP Symbolic Logic ---
# ToDo: Work with cache?

function fp_symbolic_expr(next_moment_degree)
    @variables theta phi
    
    # 1. Full Rotation (Degrees 0 to N)
    # We still need the full matrix/vector to calculate the rotation correctly
    rot_matrix = rot(next_moment_degree, theta, phi)
    vars_full  = vars_symbolic(next_moment_degree)
    
    rot_mom_list = rot_matrix * vars_full
    
    # 2. Select the specific moment expression
    indices = nidx(next_moment_degree)
    last_rot_mom = rot_mom_list[indices[end]] # The head of the shell N
    
    # 3. Extract Coefficients
    # FIX: Use ONLY the variables for the current shell (N=5 gives 4 vars)
    unknowns_shell = get_unknowns_shell(next_moment_degree)
    
    coeffs = Symbolics.gradient(last_rot_mom, unknowns_shell)
    
    return coeffs, theta, phi
end
fp_symbolic_expr(2) # Test call

# REVISED: CompileFP
function compile_fp(degree::Int)
    coeffs_sym, th, ph = fp_symbolic_expr(degree)
    
    # Force conversion to generic expressions to avoid Type/Array errors
    expr_vector = [Symbolics.toexpr(simplify(c)) for c in coeffs_sym]
    
    # Manually construct the AST for a standard Julia Vector output
    func_code = :(( $(Symbolics.toexpr(th)), $(Symbolics.toexpr(ph)) ) -> [$(expr_vector...)])
    
    func = eval(func_code)
    
    return func
end

# Wrapper to match Mathematica's FPOptimized immediately
function fp_optimized(degree, theta_val, phi_val)
    f = compile_fp(degree)
    # We now call it knowing it is the standard function f(theta, phi)
    return f(theta_val, phi_val)
end

# Just some wrapper for HyQMOM closure call
function run_gramian_closure(mom_list)
    # 1. Setup Configuration
    M = length(mom_list)-1 # Or specific logic for M
    Kn = 1.0
    closure_type = "ExtGram"
    
    # 2. Define Equations (Using your class/struct)
    equations = GramianMomentEquations1D(M, Kn, closure_type)
    
    # 3. Call the external library
    # 'u' in your snippet corresponds to 'mom_list'
    trixi_closure = HyQMOM.closure(mom_list, equations)
    
    return trixi_closure
end



########################################
################ For IC ################
########################################
# --- 2. Physics & Moment Calculations ---
# Parameters for the Maxwellian Mixture
# Important for IC
const PARAMS = (
    rho1 = 0.4,
    v1 = 0.0,
    theta1 = 0.6,
    rho2 = 0.6,
    theta2 = 0.6,
    theta3 = 0.6
)

# Helper: Calculate E[Z^k] for Standard Normal Z ~ N(0,1)
function std_normal_moment(k::Int)
    if isodd(k)
        return 0.0
    else
        # E[Z^k] = (k-1)!! which is equivalent to factorial(k) / (2^(k/2) * factorial(k/2))
        return factorial(k) / (2^(k/2) * factorial(div(k, 2)))
    end
end

# Helper: Calculate E[X^n] for X ~ N(mu, sigma)
function raw_moment_normal(mu, sigma, n)
    # Expansion of E[(mu + sigma*Z)^n]
    total = 0.0
    for k in 0:n
        # Binomial coefficient: binomial(n, k)
        coef = binomial(n, k)
        term = coef * (mu^(n - k)) * (sigma^k) * std_normal_moment(k)
        total += term
    end
    return total
end

# REVISED: Moment1D
function moment_1d(mu, sigma2, n)
    # We pass standard deviation (sqrt of variance) to our helper
    return raw_moment_normal(mu, sqrt(sigma2), n)
end

# REVISED: Moment3D (No changes needed here, just calling the fixed moment_1d)
function moment_3d(indices, v2)
    i, j, k = indices
    p = PARAMS
    
    term1 = p.rho1 * moment_1d(p.v1, p.theta1, i) * moment_1d(0.0, p.theta3, j) * moment_1d(0.0, p.theta3, k)
    term2 = p.rho2 * moment_1d(v2, p.theta2, i) * moment_1d(0.0, p.theta3, j) * moment_1d(0.0, p.theta3, k)
            
    return term1 + term2
end

# Equivalent to: MomList[v2_, degree_]
function mom_list(v2, degree)
    # Map the Moment3D function over the generated indices
    indices = mainmomindex(degree)
    return [moment_3d(idx, v2) for idx in indices]
end

# Equivalent to: NextOrderMomList[v2_, nextDegree_]
function next_order_mom_list(v2, next_degree)
    indices = index(next_degree)
    return [moment_3d(idx, v2) for idx in indices]
end

# --- 3. Example Execution (Optional) ---
# To see the result, you can run:
println(mom_list(1.5, 2))
println(next_order_mom_list(0.5, 4))


# --- M=4 Case Setup ---
M = 4
# 1. Compile FP for Degree 5 (as requested by FPM4 = CompileFP[5])
# We store the function, though it is not directly used in the bvec loop below.
const FPM4_FUNC = compile_fp(M+1);

# 2. Define Angles
const ANGLES_M4 = [
    (3.14159, 1.5708),
    (0.684719, 4.71239),
    (2.03444, 1.5708),
    (2.18628, 0.886077), # ! The last angle seems to make problems
    # (0.0, 0.0)
];

# Run the compiled FP function on the angles (resultM4)
result_M4 = [FPM4_FUNC(theta, phi) for (theta, phi) in ANGLES_M4];
A_matrix = hcat(result_M4...)';
display(A_matrix)
# ! A_matrix is correct

# Example run how I could use it in Trixi.jl
v2 = 1.5; # Example velocity
M = 4
moments_init = mom_list(v2, M) # Init moments of degree 4
# moments_init = [1., 0.9, 0, 0, 1.95, 0., 0., 0.6, 0, 0.6, 3.645, 0, 0, 0.54, 0., 0.54, 0, 0, 0, 0, 8.9775, 0., 0., 1.17, 0, 1.17, 0., 0., 0., 0., 1.08, 0, 0.36, 0, 1.08]
# Set everything to zero where index_1d says so
# ToDo: Something like this, but not this, not working :(
# ToDo: Find only the non-zero indices directly
valid_indices = []
index_start = 0
for i in 1:M
    append!(valid_indices, index_1d(i) .+ nidx(M)[i] .- 1) # Offset by shell start
end
target_indices = nidx(M) # What we use for the closure
rhs = Num[] # To store the rhs

for (theta, phi) in ANGLES_M4
    # 3. Get Rotation Matrix (Using NEW idx order)
    R = rot(M, theta, phi)

    # 4. Rotate
    rotated_moments_full = R * moments_init

    # 5. Extract (m0, m1, m2, m3, m4)
    substituted = rotated_moments_full[target_indices]

    # 6. Closure
    equations = GramianMomentEquations1D(M, 1.0, "ExtGram") # Kn and ExtGram currently "placeholder"
    trixi_closure = HyQMOM.closure(substituted, equations)
    
    push!(rhs, trixi_closure)
end
# Solve system Ax = b where b is 'rhs' and A is from A_matrix
# ? The last angle of rhs is calculated incorrectly, but why is this?
approx_moments = A_matrix \ rhs;

# Get Exact Result
exact_moments = get_exact_next_order_moments(v2, M+1);

println("\n--- Example Closure Result M4 (Approximated Moments) ---")
println("Predicted moments: ", approx_moments)
println("Exact moments: ", exact_moments)
println("Relative Error: ", abs.(exact_moments .- approx_moments) ./ abs.(exact_moments))




# --- 4. bvecM4 Calculation ---
# This is where the magic happens
function bvec_m4(v2)
    # 1. Calculate Raw Moments (Degree 4)
    raw_moments = mom_list(v2, M)

    # 2. Get Indices for Shell Heads (1, 2, 5, 11, 21)
    target_indices = nidx(M)

    results = Num[] # "Num" for symbolic compatibility, # ToDo: Change to Float64 later?

    for (theta, phi) in ANGLES_M4
        # 3. Get Rotation Matrix (Using NEW idx order)
        R = rot(M, theta, phi)

        # 4. Rotate
        rotated_moments_full = R * raw_moments

        # 5. Extract (m0, m1, m2, m3, m4)
        substituted = rotated_moments_full[target_indices]

        # 6. Closure
        val = run_gramian_closure(substituted)
        
        push!(results, val)
    end
    
    return results
end

# --- Execution ---

# Define v2list (Example range, as it wasn't defined in the snippet)
v2_list = [0.1, 0.3, 0.8, 1.5, 2.5, 3.0];

# Generate bM4 Table
println("Calculating bM4 for ", length(v2_list), " velocity points...")
b_M4 = [bvec_m4(v) for v in v2_list];

# Display Result (formatted like MatrixForm)
# b_M4 is a Vector of Vectors. We convert to Matrix for display.
b_M4_matrix = hcat(b_M4...)'; # Transpose to get rows=v2, cols=angles
display(b_M4_matrix)
# ! The last angle seems to make problems






#################################################
######## Closure & Error Calculation ############
#################################################

# ! Not necessary for the Trixi implementation, but kept for completeness


# --- 1. Helper: Get Exact Values for the Shell ---
# This calculates the theoretical moments and filters them exactly like 'get_unknowns_shell'
# Matches: NextOrderMomList[v2, 5][[index1D[5]]]
# Equivalent to Mathematica: NextOrderMomList[v2, 5][[index1D[5]]]
function get_exact_next_order_moments(v2, degree)
    # 1. Get indices ONLY for the specific shell (sum = degree)
    shell_indices = index(degree) 
    
    # 2. Calculate the moment value for each index in this shell
    # This ensures we get Degree 5 values, not Degree 0 or 1
    shell_values = [moment_3d(idx, v2) for idx in shell_indices]
    
    # 3. Filter to keep only the "Unknowns" (Symmetric moments)
    # Matches the logic of index1D / generate1D
    filtered_values = Float64[]
    
    for (k, (i, j, m)) in enumerate(shell_indices)
        # Check the symmetry condition: Even y, Even z, y >= z
        if iseven(j) && iseven(m) && (j >= m)
            push!(filtered_values, shell_values[k])
        end
    end
    
    return filtered_values
end

# --- 2. Setup Matrices ---

# Construct Matrix 'A' from result_M4
# result_M4 is a Vector of Vectors. 
# We stack them to create the matrix A where A * x = b
# Transpose is needed because hcat stacks as columns, but result_M4 are rows of coeffs
A_matrix = hcat(result_M4...)';

println("Linear System Dimensions:")
println("Matrix A: ", size(A_matrix))

# --- 3. Main Loop: Solve and Error Calculation ---

closure_result_M4 = Vector{Num}[]; # "Num" for symbolic compatibility, # ToDo: Change to Float64 later?
exact_result_M4 = Vector{Float64}[];
RE_M4 = Vector{Num}[]; # "Num" for symbolic compatibility, # ToDo: Change to Float64 later?

# b_M4_cols contains the b-vectors for each v2 (from previous step)
for (i, v2) in enumerate(v2_list)
    # 1. Get RHS vector b
    b_vec = b_M4[i]
    
    # 2. Solve Ax = b
    # The backslash operator (\) is Julia's efficient linear solver
    x_approx = A_matrix \ b_vec
    push!(closure_result_M4, x_approx)
    
    # 3. Get Exact Result
    x_exact = get_exact_next_order_moments(v2, 5)
    push!(exact_result_M4, x_exact)
    
    # 4. Calculate Relative Error
    # Abs[Exact - Closure] / Abs[Exact]
    # We use a small epsilon to avoid division by zero if Exact is 0
    re_vec = abs.(x_exact .- x_approx) ./ (abs.(x_exact) .+ 1e-16)
    push!(RE_M4, re_vec)
end

# --- 4. RMS Calculation ---
# RERMSM4 = Table[Sqrt[Mean[(REM4[[i]])^2]], {i, 1, 6}]
RERMS_M4 = [sqrt(mean(re.^2)) for re in RE_M4]

# --- 5. Display Results ---

println("\n--- Closure Result M4 (Approximated Moments) ---")
# Convert Vector of Vectors to Matrix for display
display(hcat(closure_result_M4...)')

println("\n--- Exact Result M4 ---")
display(hcat(exact_result_M4...)')

println("\n--- Relative Error (RMS) per Velocity Step ---")
# Display as column vector
display(reshape(RERMS_M4, :, 1))