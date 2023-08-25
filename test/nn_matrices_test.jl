module nn_matrices_test

proj_dir =  dirname(dirname(@__FILE__))*"\\src\\"
include(proj_dir*"nn_matrices.jl")
using Test

@testset "make_NN_mat" begin
    for i ∈ 1:5
        m,n = rand(1:150,2)
        M = nn_matrices.make_NN_mat(m,n)
        @test size(M) == (m,n)
        @test all(M .≥ 0)
    end
end


@testset "testNN" begin
    M = rand(10,10)
    @test testNN(M) == true
end


end
