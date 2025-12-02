# Disclaimer: This Julia code is a translation of a Mathematica script 
# The Translation was mainly done with Google Gemini and manually refined.
# Some constructs may not have direct equivalents and were adapted accordingly.

using Symbolics
using LinearAlgebra
using Distributions
using Latexify
if !endswith(Base.active_project(), "../Project.toml")
    import Pkg; Pkg.activate(".")
end # Runs in environment setup
using HyQMOM, Trixi

# --- 1. Index Generation Functions ---

# Equivalent to: idx[n_]
function idx(n)
    # 1. Generate tuples
    # 2. collect(...) creates a Multidimensional Array
    # 3. vec(...) flattens it into a 1D Vector, which unique requires
    raw_tuples = vec(collect(Iterators.product(fill(1:3, n)...)))
    
    # Convert tuples to vectors and sort the elements inside each vector
    # e.g., (3, 1, 2) -> [1, 2, 3]
    sorted_vectors = [sort(collect(t)) for t in raw_tuples]
    
    # Remove duplicates
    return unique(sorted_vectors)
end

# Equivalent to: index[n_]
# ! needed
function index(n)
    # Generate tuples 0 to n of length 3
    # Iterators.product creates the Cartesian product
    raw_tuples = collect(Iterators.product(0:n, 0:n, 0:n))
    
    # Flatten and convert to vectors
    vecs = [collect(t) for t in raw_tuples]
    
    # Select those where sum is n
    selected = filter(v -> sum(v) == n, vecs)
    
    # Mathematica sorts lexicographically then Reverses. 
    # We sort standardly then reverse to match the script perfectly.
    sort!(selected)
    return reverse(selected)
end

function tuple_to_powers(tup, n)
    c1 = count(==(1), tup)
    c2 = count(==(2), tup)
    c3 = count(==(3), tup)
    return [c1, c2, c3]
end

# Equivalent to: index1D[n_]
function index_1d(n)
    # Get the standard index list
    full_indices = index(n)
    
    # Find positions where y and z are even, and y >= z
    # Return linear indices (Flatten[Position[...]])
    pos = Int[]
    for (i, val) in enumerate(full_indices)
        _, y, z = val
        if iseven(y) && iseven(z) && (y >= z)
            push!(pos, i)
        end
    end
    return pos
end

# Equivalent to: mainmomindex[nmax_]
# ! needed
function mainmomindex(nmax)
    # Flatten[Table[...], 1] is equivalent to vcat in Julia
    # We map the index function over the range 0 to nmax
    return vcat([index(n) for n in 0:nmax]...)
end

# --- 2. Physics & Moment Calculations ---

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
# ! needed
function moment_1d(mu, sigma2, n)
    # We pass standard deviation (sqrt of variance) to our helper
    return raw_moment_normal(mu, sqrt(sigma2), n)
end

# REVISED: Moment3D (No changes needed here, just calling the fixed moment_1d)
# ! needed
function moment_3d(indices, v2)
    i, j, k = indices
    p = PARAMS
    
    term1 = p.rho1 * moment_1d(p.v1, p.theta1, i) * moment_1d(0.0, p.theta3, j) * moment_1d(0.0, p.theta3, k)
    term2 = p.rho2 * moment_1d(v2, p.theta2, i) * moment_1d(0.0, p.theta3, j) * moment_1d(0.0, p.theta3, k)
            
    return term1 + term2
end

# Equivalent to: MomList[v2_, degree_]
# ! needed
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



# --- 3. Tensor Transformation Logic ---

# Memoization Cache: Key is (degree, theta_symbol, phi_symbol)
const SYMBOLIC_CACHE = Dict{Any, Matrix{Num}}()

function get_rotation_matrix(theta, phi)
    # Using standard Julia cos/sin. 
    # If theta/phi are symbolic, these return symbolic expressions.
    ct, st = cos(theta), sin(theta)
    cp, sp = cos(phi), sin(phi)
    
    return [
        -ct*sp  -st*sp   cp;
         st*sp  -ct*sp   0;
        -ct*cp  -cp*st  -sp
    ]
end

function tensor_transformation_symbolic(n, theta, phi)
    # Check cache
    key = (n, theta, phi)
    if haskey(SYMBOLIC_CACHE, key)
        return SYMBOLIC_CACHE[key]
    end

    row_basis = idx(n) 
    col_basis = index(n)
    col_lookup = Dict(vec => i for (i, vec) in enumerate(col_basis))
    
    # CRITICAL CHANGE: Initialize with symbolic zeros (Num(0))
    # We use size matching row/col basis
    dim_rows = length(row_basis)
    dim_cols = length(col_basis)
    TT = zeros(Num, dim_rows, dim_cols)
    
    R = get_rotation_matrix(theta, phi)
    all_tuples = collect(Iterators.product(fill(1:3, n)...))

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

    # Simplify the resulting matrix to group sin/cos terms
    # broadcasting simplify(.) over the matrix
    TT_simplified = simplify.(TT)
    
    SYMBOLIC_CACHE[key] = TT_simplified
    return TT_simplified
end

# 1. Create a specific cache for numeric matrices to speed up loops
# ! needed
const NUMERIC_CACHE = Dict{Tuple{Int, Float64, Float64}, Matrix{Float64}}()

# ! needed
function tensor_transformation_memo(n::Int, theta::Real, phi::Real)
    # Ensure inputs are Float64 for the cache key
    val_theta = Float64(theta)
    val_phi = Float64(phi)
    key = (n, val_theta, val_phi)
    
    # 1. Check Numeric Cache
    if haskey(NUMERIC_CACHE, key)
        return NUMERIC_CACHE[key]
    end

    # 2. Call the Symbolic Function
    # Even if we pass numbers, this returns Matrix{Num} because of initialization
    mat_num_wrapped = tensor_transformation_symbolic(n, val_theta, val_phi)
    
    # 3. Unwrap and Convert to Float64
    # Symbolics.value unwraps the 'Num' type to the underlying number
    mat_float = Float64.(Symbolics.value.(mat_num_wrapped))
    
    # 4. Save and Return
    NUMERIC_CACHE[key] = mat_float
    return mat_float
end


# --- 3. Run Comparison ---

# Define symbolic variables
@variables θ φ 

# Calculate the transformation matrix for n=2
TT_sym = tensor_transformation_symbolic(2, θ, φ)
latexify(TT_sym)

# Display result
# This prints the symbolic matrix to the terminal
display(TT_sym)




# --- 4. Symbolic Block Diagonal Construction ---

# Cache for the final block matrices
const ROT_CACHE = Dict{Any, Matrix{Num}}()

# Helper function to manually construct a Block Diagonal Matrix
# (Stitches a list of matrices [M1, M2...] into one large matrix)
function construct_block_diagonal(matrices::Vector{Matrix{Num}})
    # Calculate total dimensions
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)
    
    # Initialize large symbolic matrix with zeros
    # We use 'Num' to ensure it can hold symbolic variables
    full_matrix = zeros(Num, total_rows, total_cols)
    
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
function rot_symbolic(moment_degree, theta, phi)
    # Check cache
    key = (moment_degree, theta, phi)
    if haskey(ROT_CACHE, key)
        return ROT_CACHE[key]
    end

    # Collect the blocks for n = 0 to moment_degree
    blocks = Matrix{Num}[]
    
    for n in 0:moment_degree
        if n == 0
            # IdentityMatrix[1] is just a 1x1 matrix containing 1
            # We explicitly make it 'Num' type to match the rest
            push!(blocks, reshape([Num(1)], 1, 1))
        else
            # Call the tensor transformation function we defined earlier
            mat = tensor_transformation_symbolic(n, theta, phi)
            push!(blocks, mat)
        end
    end
    
    # Stitch them together
    final_matrix = construct_block_diagonal(blocks)
    
    # Save and return
    ROT_CACHE[key] = final_matrix
    return final_matrix
end

# --- 5. Final Execution ---

@variables θ φ 

# Calculate the full block diagonal matrix for Degree 2
# This combines n=0 (scalar), n=1 (3x3), and n=2 (6x6)
M = rot_symbolic(2, θ, φ)

println("Generated Block Diagonal Matrix of size: ", size(M))
display(M)



# --- Helper: Numeric Block Diagonal Construction ---
# ! needed
function construct_numeric_block_diagonal(matrices::Vector{Matrix{Float64}})
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)
    
    full_matrix = zeros(Float64, total_rows, total_cols)
    
    row_cursor = 1
    col_cursor = 1
    
    for m in matrices
        r, c = size(m)
        full_matrix[row_cursor:(row_cursor + r - 1), 
                    col_cursor:(col_cursor + c - 1)] = m
        row_cursor += r
        col_cursor += c
    end
    
    return full_matrix
end

# --- 1. RotNumeric Implementation ---
# Equivalent to: RotNumeric[MomentDegree_, theta, phi]
# ! needed
function rot_numeric(moment_degree, theta, phi)
    blocks = Matrix{Float64}[]
    
    for n in 0:moment_degree
        if n == 0
            # IdentityMatrix[1] as Float64
            push!(blocks, reshape([1.0], 1, 1))
        else
            # Call the NUMERIC tensor transformation (from previous steps)
            # We assume tensor_transformation_memo (or tensor_transformation) exists
            mat = tensor_transformation_memo(n, Float64(theta), Float64(phi))
            push!(blocks, mat)
        end
    end
    
    return construct_numeric_block_diagonal(blocks)
end

# --- 2. VarsSymbolic Implementation ---
# Equivalent to: VarsSymbolic[MomentDegree]
# This generates the vector of variables V_i_j_k, filtering for the 1D symmetries
function vars_symbolic(moment_degree)
    # Initialize an empty vector of Symbolic Numbers
    vars = Num[]
    
    for n in 0:moment_degree
        if n == 0
            # n=0 is always valid (V_0_0_0)
            push!(vars, first(@variables V_0_0_0))
        else
            # 1. Get the FULL list of indices for this degree
            # (Matches the length of the rotation matrix columns)
            full_indices = index(n) 
            
            for (i, j, k) in full_indices
                # 2. Check the validity condition:
                #    EvenQ[y] && EvenQ[z] && (y >= z)
                #    Here: y->j, z->k
                if iseven(j) && iseven(k) && (j >= k)
                    # Case A: Valid -> Generate Symbol
                    # We name it V_i_j_k to be Julia-safe
                    sym_name = Symbol("V_$(i)_$(j)_$(k)")
                    var = Symbolics.unwrap(first(@variables $sym_name))
                    push!(vars, Num(var))
                else
                    # Case B: Invalid -> Insert Zero
                    push!(vars, Num(0))
                end
            end
        end
    end
    
    return vars
end

# --- Execution ---

# Generate the vector
V_sym = vars_symbolic(2)

# Display as a column vector (MatrixForm equivalent)
display(reshape(V_sym, :, 1))



# Equivalent to: nIDX[MomentDegree_] := Table[1 + 1/6 k (k + 1) (k + 2), ...]
# ! needed
function nidx(moment_degree)
    # Julia is 1-based, matching the specific logic of the Mathematica script
    # This generates the start index of every moment shell in the flattened list
    return Int[1 + div(k*(k+1)*(k+2), 6) for k in 0:moment_degree]
end

# --- 2. Unknowns Helper (x) ---
# Equivalent to: x[NextMomentdegree_] := DeleteCases[generate1D[NextMomentdegree], 0]
# Equivalent to Mathematica's: DeleteCases[generate1D[n], 0]
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

# --- 3. FP Symbolic Logic ---

# Cache for the compiled functions
# ! needed
const FP_CACHE = Dict{Int, Function}()

# Equivalent to: FPSymbolic[NextMomentdegree_]
# ! needed
function fp_symbolic_expr(next_moment_degree)
    @variables theta phi
    
    # 1. Full Rotation (Degrees 0 to N)
    # We still need the full matrix/vector to calculate the rotation correctly
    rot_matrix = rot_symbolic(next_moment_degree, theta, phi)
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

# REVISED: CompileFP
# ! needed
function compile_fp(degree::Int)
    if haskey(FP_CACHE, degree)
        return FP_CACHE[degree]
    end

    coeffs_sym, th, ph = fp_symbolic_expr(degree)
    
    # Force conversion to generic expressions to avoid Type/Array errors
    expr_vector = [Symbolics.toexpr(simplify(c)) for c in coeffs_sym]
    
    # Manually construct the AST for a standard Julia Vector output
    func_code = :(( $(Symbolics.toexpr(th)), $(Symbolics.toexpr(ph)) ) -> [$(expr_vector...)])
    
    func = eval(func_code)
    
    FP_CACHE[degree] = func
    return func
end

# Wrapper to match Mathematica's FPOptimized immediately
function fp_optimized(degree, theta_val, phi_val)
    f = compile_fp(degree)
    # We now call it knowing it is the standard function f(theta, phi)
    return f(theta_val, phi_val)
end

# The function wrapper you requested
# ! needed
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

# --- 5. Execution Test ---

println("--- Testing nIDX ---")
println("Indices for Degree 2: ", nidx(2))

println("\n--- Testing Unknowns (x) ---")
display(get_unknowns_shell(2))

# Clean cache if needed (optional)
empty!(FP_CACHE)

println("\n--- Testing FPSymbolic (Optimization) ---")
fast_fp_func = compile_fp(2)

# Run it
theta_val = 0.5
phi_val = 0.3
result_coeffs = fast_fp_func(theta_val, phi_val)

println("Success! Coefficients:")
display(result_coeffs)

println("\n--- Testing HyQMOM Wrapper ---")
# Dummy list for testing
dummy_moments = rand(5) 
run_gramian_closure(dummy_moments)



# --- M=4 Case Setup ---

# 1. Compile FP for Degree 5 (as requested by FPM4 = CompileFP[5])
# We store the function, though it is not directly used in the bvec loop below.
const FPM4_FUNC = compile_fp(5)

# 2. Define Angles
const ANGLES_M4 = [
    (3.14159, 1.5708),
    (0.684719, 4.71239),
    (2.03444, 1.5708),
    (2.18628, 0.886077), # ! The last angle seems to make problems
]

# Run the compiled FP function on the angles (resultM4)
result_M4 = [FPM4_FUNC(theta, phi) for (theta, phi) in ANGLES_M4]
display(result_M4)

# 3. Define Constants
const NEXT_M4 = 5
const MD_M4 = NEXT_M4 - 1 # Degree 4

# 4. Define bvecM4 Function
# Logic: Rotate moments numerically, pick specific indices, run closure.
function bvec_m4(v2)
    # Constants
    md_m4 = 4 # Degree 4
    
    # 1. Calculate Raw Moments (Degree 4)
    raw_moments = mom_list(v2, md_m4)

    # 2. Get Indices for Shell Heads (1, 2, 5, 11, 21)
    target_indices = nidx(md_m4)

    results = Float64[]

    for (theta, phi) in ANGLES_M4
        # 3. Get Rotation Matrix (Using NEW idx order)
        R = rot_numeric(md_m4, theta, phi)

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
# v2_list = collect(0.0:0.2:1.0) 
v2_list = [0.1, 0.3, 0.8, 1.5, 2.5, 3.0]

# Generate bM4 Table
println("Calculating bM4 for ", length(v2_list), " velocity points...")
b_M4 = [bvec_m4(v) for v in v2_list]

# Display Result (formatted like MatrixForm)
# b_M4 is a Vector of Vectors. We convert to Matrix for display.
b_M4_matrix = hcat(b_M4...)' # Transpose to get rows=v2, cols=angles
display(b_M4_matrix)
# ! The last angle seems to make problems




# --- 1. Helper: Get Exact Values for the Shell ---
# This calculates the theoretical moments and filters them exactly like 'get_unknowns_shell'
# Matches: NextOrderMomList[v2, 5][[index1D[5]]]
# Equivalent to Mathematica: NextOrderMomList[v2, 5][[index1D[5]]]
# ! needed for comparison
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
A_matrix = hcat(result_M4...)' 

println("Linear System Dimensions:")
println("Matrix A: ", size(A_matrix))

# --- 3. Main Loop: Solve and Error Calculation ---

closure_result_M4 = Vector{Float64}[]
exact_result_M4 = Vector{Float64}[]
RE_M4 = Vector{Float64}[]

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