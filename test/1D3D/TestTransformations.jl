idx4_true = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 3],
    [1, 1, 2, 2],
    [1, 1, 2, 3],
    [1, 1, 3, 3],
    [1, 2, 2, 2],
    [1, 2, 2, 3],
    [1, 2, 3, 3],
    [1, 3, 3, 3],
    [2, 2, 2, 2],
    [2, 2, 2, 3],
    [2, 2, 3, 3],
    [2, 3, 3, 3],
    [3, 3, 3, 3]
]
@test HyQMOM.idx(4) == idx4_true

index4_true = [
    [4, 0, 0],
    [3, 1, 0],
    [3, 0, 1],
    [2, 2, 0],
    [2, 1, 1],
    [2, 0, 2],
    [1, 3, 0],
    [1, 2, 1],
    [1, 1, 2],
    [1, 0, 3],
    [0, 4, 0],
    [0, 3, 1],
    [0, 2, 2],
    [0, 1, 3],
    [0, 0, 4]
]

@test HyQMOM.index(4) == index4_true

index_1d_4_true = [1, 4, 11, 13]
@test HyQMOM.index_1d(4) == index_1d_4_true

mainmomindex4_true = [
    [0, 0, 0], # U
    [1, 0, 0], # U_x
    [0, 1, 0],
    [0, 0, 1],
    [2, 0, 0], # U_xx
    [1, 1, 0],
    [1, 0, 1],
    [0, 2, 0],
    [0, 1, 1],
    [0, 0, 2],
    [3, 0, 0], # U_xxx
    [2, 1, 0],
    [2, 0, 1],
    [1, 2, 0],
    [1, 1, 1],
    [1, 0, 2],
    [0, 3, 0],
    [0, 2, 1],
    [0, 1, 2],
    [0, 0, 3],
    [4, 0, 0], # U_xxxx
    [3, 1, 0],
    [3, 0, 1],
    [2, 2, 0],
    [2, 1, 1],
    [2, 0, 2],
    [1, 3, 0],
    [1, 2, 1],
    [1, 1, 2],
    [1, 0, 3],
    [0, 4, 0],
    [0, 3, 1],
    [0, 2, 2],
    [0, 1, 3],
    [0, 0, 4]
]
@test HyQMOM.mainmomindex(4) == mainmomindex4_true

# Symbolic variables
rotation_matrix = [
    [-0.259343, -0.14168, 0.955336],
    [0.14168, -0.259343, 0.0],
    [-0.838387, -0.458013, -0.29552]
]
# @test HyQMOM.get_rotation_matrix(0.5, 0.3) .≈ rotation_matrix atol=1e-5

tensor_transformation4_true = [
    0.067259 0.0734875 -0.49552 0.0200732 -0.270704 0.912668;
    -0.0367438 0.0471858 0.135352 0.0367438 -0.24776 0.0;
    0.21743 0.237565 -0.7243 0.0648912 -0.395687 -0.282321;
    0.0200732 -0.0734875 0.0 0.067259 0.0 0.0;
    -0.118783 0.152539 -0.0418693 0.118783 0.0766412 0.0;
    0.702892 0.767983 0.49552 0.209776 0.270704 0.0873322
]

@test isapprox(HyQMOM.tensor_transformation(2, 0.5, 0.3), tensor_transformation4_true; atol=1e-5, rtol=1e-5)

rotation2_true = [
    1 0 0 0 0 0 0 0 0 0;
    0 -0.259343 -0.14168 0.955336 0 0 0 0 0 0;
    0 0.14168 -0.259343 0.0 0 0 0 0 0 0;
    0 -0.838387 -0.458013 -0.29552 0 0 0 0 0 0;
    0 0 0 0 0.067259 0.0734875 -0.49552 0.0200732 -0.270704 0.912668;
    0 0 0 0 -0.0367438 0.0471858 0.135352 0.0367438 -0.24776 0.0;
    0 0 0 0 0.21743 0.237565 -0.7243 0.0648912 -0.395687 -0.282321;
    0 0 0 0 0.0200732 -0.0734875 0.0 0.067259 0.0 0.0;
    0 0 0 0 -0.118783 0.152539 -0.0418693 0.118783 0.0766412 0.0;
    0 0 0 0 0.702892 0.767983 0.49552 0.209776 0.270704 0.0873322
]

@test isapprox(HyQMOM.rot(2, 0.5, 0.3), rotation2_true; atol=1e-5, rtol=1e-5)


nidx4_true = [1, 2, 5, 11, 21]
@test HyQMOM.nidx(4) == nidx4_true

# Test case
M = 4
FPM4_FUNC = HyQMOM.compile_fp(M+1);

ANGLES_M4 = [
    (3.14159, 1.5708),
    (0.684719, 4.71239),
    (2.03444, 1.5708),
    (2.18628, 0.886077)
];

# Run the compiled FP function on the angles (resultM4)
result_M4 = [FPM4_FUNC(theta, phi) for (theta, phi) in ANGLES_M4];
A_matrix = hcat(result_M4...)';

A_matrix_true = [
    1.0 7.04154e-11 2.47916e-22 2.85023e-21;
    0.278855 1.85903 0.619677 9.66339e-12;
    0.0178878 0.715528 1.43108 1.44815e-10;
    0.017889 0.357775 0.357769 2.14663
]

@test isapprox(A_matrix, A_matrix_true; atol=1e-5, rtol=1e-5)

let # Start a local scope
    v2 = 1.5; # Example velocity
    moments_init = [1., 0.9, 0, 0, 1.95, 0., 0., 0.6, 0, 0.6, 3.645, 0, 0, 0.54, 0., 0.54, 0, 0, 0, 0, 8.9775, 0., 0., 1.17, 0, 1.17, 0., 0., 0., 0., 1.08, 0, 0.36, 0, 1.08]#HyQMOM.mom_list(v2, M) # Init moments of degree 4
    slab_indices = []
    index_start = 0
    for i in 0:M
        append!(slab_indices, HyQMOM.index_1d(i) .+ index_start) # Offset by shell start
        index_start += size(HyQMOM.index(i), 1)
    end
    moments_init = [i in slab_indices ? moments_init[i] : 0.0 for i in 1:length(moments_init)]
    target_indices = HyQMOM.nidx(M) # What we use for the closure
    using Symbolics # for Num
    # ToDo: Only use Float64
    rhs = Num[] # To store the rhs

    Kn = 1.0
    closure_type = "ExtGram"
    equations = GramianMomentEquations1D3D(M, Kn, closure_type)

    for (theta, phi) in ANGLES_M4
        # 3. Get Rotation Matrix (Using NEW idx order)
        R = HyQMOM.rot(M, theta, phi)

        # 4. Rotate
        rotated_moments_full = R * moments_init

        # 5. Extract (m0, m1, m2, m3, m4)
        substituted = rotated_moments_full[target_indices]

        val = HyQMOM.closure(substituted, equations)
            
        push!(rhs, val)
    end
    approx_moments = A_matrix \ rhs;

    rhs4_true = [21.62986080236905, 10.739459101432965, 3.3612664962508125, 2.1945213576183544]

    @test isapprox(rhs, rhs4_true; atol=1e-5, rtol=1e-5)

    approx_moments4_true = [21.62986080332397, 2.207547005105253, 0.974641967316188, 0.31169070027773604]

    @test isapprox(approx_moments, approx_moments4_true; atol=1e-5, rtol=1e-5)

    # Test the closure_transform implementation
    moments_init_reduced = moments_init[slab_indices]
    @test isapprox(approx_moments4_true, HyQMOM.closure_transform(moments_init_reduced, equations); atol=1e-5, rtol=1e-5)
end