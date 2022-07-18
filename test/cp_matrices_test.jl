module moments_test
using Test
using LinearAlgebra

src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"cp_matrices.jl")
include(src_dir*"matrix_IO.jl")
const cpm = cp_matrices

@testset "generate_lit_cp_mats" begin
  lit_cp_mats = cpm.generate_lit_cp_mats()
  K = [keys(lit_cp_mats)...]
  for k in K
    M = lit_cp_mats[k]
    @test cpm.isPSD(M) 
    @test cpm.isNN(M)  
    @test cpm.isDNN(M) 
    @test cpm.isDiagOne(M) 

    n, nzd, r, _ = matrix_IO.split_name(k) 
    @test size(M)[1] == n
    @test rank(M) == r
    @test isapprox(nzd, round((sum(M .> 0) - n)/2/binomial(round(Int,n),2),digits=2), atol=0.02)
  end
end

@testset "generate_lit_non_cp_mats" begin
  mats = cp_matrices.generate_lit_non_cp_mats()
  K = [keys(mats)...]
  for k in K
    M = mats[k]
    @test cpm.isPSD(M) 
    @test cpm.isNN(M)  
    @test cpm.isDNN(M) 
    @test cpm.isDiagOne(M) 

    n, nzd, r, _ = matrix_IO.split_name(k) 
    @test size(M)[1] == n
    @test rank(M) == r
    @test isapprox(nzd, round((sum(M .> 0) - n)/2/binomial(round(Int,n),2),digits=2), atol=0.02)
  end
end

@testset "generate_random_cp_mats" begin
  N = [5,6,7,8,9]
  L = vcat(map(k->div(binomial(k,2),2),[5,6,7,8,9])...)
  K = ones(Int,length(L))*3
  mats, rand_cp_facts = cpm.generate_random_cp_mats(N, L, K)
  K = [keys(mats)...]
  for k in K
    println(k)
    M = mats[k]
    @test cpm.isPSD(M) 
    @test cpm.isNN(M)  
    @test cpm.isDNN(M) 
    @test cpm.isDiagOne(M) 

    n, nzd, r, ur = matrix_IO.split_name(k) 
    @test size(M)[1] == n
    @test rank(M) == r
    @test length(rand_cp_facts["R_"*k]) == ur
    @test isapprox(nzd, round((sum(M .> 0) - n)/2/binomial(round(Int,n),2),digits=2), atol=0.02)
  end
end


end

# function CP_necessities(A) 
#   @test isPSD(A)
#   @test isNN(A) 
#   @test isDiagOne(A)
# end 
  
# function run_tests()
#   @testset "Misc." begin
#     @test circmatrix([1,2,3])  == [3 2 1; 1 3 2; 2 1 3]
#     @test makediagone([1,2,3,4]*[1 2 3 4]) == [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
#   end

#   @testset "get_lit_cp_mats" begin 
#     lit_cp_mats = get_lit_cp_mats()  
#     @test length(lit_cp_mats) == 27
#     for k in keys(lit_cp_mats)
#         CP_necessities(lit_cp_mats[k]) 
#     end
#   end

#   @testset "load_random_cp_mats" begin
#     load_dir = dirname(dirname(@__FILE__))*"\\assets\\data\\randomly_generated_cp-matrices\\mats\\"
#     random_cp_mats = load_random_cp_mats(load_dir)
#     @test length(random_cp_mats) == 50
#     for k in keys(random_cp_mats)
#       CP_necessities(random_cp_mats[k]) 
#     end
#   end

# end

  # @testset "get_random_cp_mats(n,r)" begin
  #     for n in 5:10, r in rand(1:20,5)
  #         r,R,A = get_random_cp_mats(n,r)
  #         @test isNN(A) 
  #         @test size(A)[1] == n
  #         @test isPSD(A)
  #     end
  # end  

  
# m=3
# n=4
# l=2 
# k=4

# @assert sum(cp_matrices.gen_random_supp_mat(m,n,l)) == l
# @assert Set(unique(cp_matrices.gen_random_supp_mat(m,n,l))) == Set([0 1])
# @assert sum(cp_matrices.gen_random_supp_mat(n,l)) == 2l+n
# @assert Set(unique(cp_matrices.gen_random_supp_mat(n,l))) == Set([0 1])
# M, fact = cp_matrices.get_random_cp_mat((n,l),k) 
# @assert length(fact) == k
# @assert cp_matrices.makediagone(sum([f*f' for f in fact])) == M

# mats_facts = cp_matrices.get_random_cp_mats(n,k)
# @assert length(mats_facts) == binomial(k,2)
