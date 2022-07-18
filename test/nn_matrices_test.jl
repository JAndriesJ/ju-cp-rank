module nn_matrices_test

proj_dir =  dirname(dirname(@__FILE__))*"\\src\\"
include(proj_dir*"nn_matrices.jl")
using Test

@testset "make_NN_mat" begin
    for i ∈ 1:5
    m,n = rand(1:150,2)
    M = nn_matrices.make_NN_mat(m*n,m,n)
    @test size(M) == (m,n)
    @test all(M .≥ 0)
    end
end

@testset "make_Gₘ_adj_mat" begin
    for i ∈ 1:5
        m,n = rand(1:150,2)
        M = nn_matrices.make_NN_mat(m*n,m,n)
        sqrtMₘₐₓ = sqrt(maximum(M))
        Gₘ_adj_mat = nn_matrices.make_Gₘ_adj_mat(M)
        @test size(Gₘ_adj_mat) == (m+n,m+n)
        @test unique(Gₘ_adj_mat) == [1.0, 0.0]
        @test Gₘ_adj_mat[1:m,1:m] == ones(m,m)
        @test Gₘ_adj_mat[m+1:m+n,m+1:m+n] == ones(n,n)
    end
end

end
