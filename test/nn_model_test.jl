src_dir =  dirname(dirname(@__FILE__))*"\\src\\"

include(src_dir*"nn_matrices.jl")
include(src_dir*"moments.jl")
include(src_dir*"nn_model.jl")
using Test
using JuMP # The solver that we use.
using MosekTools

using .moments ; const mom = moments

@testset "simple t = 1" begin
    t = 1
    M = [1 1 0 0
        1 0 1 0
        0 1 0 1
        0 0 1 1 ] .+ 0.0

    ξₜⁿⁿⁱᵈ = nn_model.modelξₜⁿⁿⁱᵈ(M,t)
    res = nn_model.computeξₜⁿⁿ(ξₜⁿⁿⁱᵈ)
    @test round(objective_value(res), digits=2) == 2.91
end

@testset "simple t = 2" begin
    t = 2
    M = [1 1 0 0
        1 0 1 0
        0 1 0 1
        0 0 1 1 ] .+ 0.0

    m,n = size(M)
    model = Model()
    @variable(model, Y[moments.make_mon_expo(m+n,2*t)] )
    Lxx = nn_model.get_Lxx(m,n,t,Y) # L([x]ₜ[x]ₜᵀ) ⪰ 0

    sec_ord_cons = nn_model.get_sec_ord_cons(M,Y) # L(xᵢxⱼ) = Mᵢⱼ
    @test length(sec_ord_cons) == sum(M .> 0)
    loc1_cons =  nn_model.get_loc1_cons(M,t,Y) # L((√Mₘₐₓxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [m + n]
    @test length(loc1_cons) == sum(m+n)
    # size(loc1_cons[7])
    #more tests needed
    loc2_cons = nn_model.get_loc2_cons(M,t,Y) # L((Mᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j ∈ Eₘ 
    @test length(loc2_cons) == sum(M .> 0)
    # more tests needed
    ideal_cons = nn_model.get_ideal_cons(M,t,Y) # L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ E s.t. Mᵢⱼ = 0 
    Gₘ_adj_mat = moments.make_Gₘ_adj_mat(M)
    mom₂ₜ₋₂ = moments.make_mon_expo(2t-2, Gₘ_adj_mat)
    @test length(ideal_cons) == length(mom₂ₜ₋₂)*sum(M .== 0.0)

    ξₜⁿⁿⁱᵈ = nn_model.modelξₜⁿⁿⁱᵈ(M,t)
    res = nn_model.computeξₜⁿⁿ(ξₜⁿⁿⁱᵈ)
    @test round(objective_value(res), digits=2) == 4.0
end

@testset "random" begin
    t = rand(1:3)
    m = rand(3:9) ; n = rand(3:5)
    M = nn_matrices.make_NN_mat(m,n)
    model = Model()
    @variable(model, Y[moments.make_mon_expo(m+n,2*t)] )
    Lxx = nn_model.get_Lxx(M,t,Y) # L([x]ₜ[x]ₜᵀ) ⪰ 0

    sec_ord_cons = nn_model.get_sec_ord_cons(M,Y) # L(xᵢxⱼ) = Mᵢⱼ
    @test length(sec_ord_cons) == sum(M .> 0)
    loc1_cons =  nn_model.get_loc1_cons(M,t,Y) # L((√Mₘₐₓxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [m + n]
    @test length(loc1_cons) == sum(m+n)
    # size(loc1_cons[7])
    #more tests needed
    loc2_cons = nn_model.get_loc2_cons(M,t,Y) # L((Mᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j ∈ Eₘ 
    @test length(loc2_cons) == sum(M .> 0)
    # more tests needed
    ideal_cons = nn_model.get_ideal_cons(M,t,Y) # L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ E s.t. Mᵢⱼ = 0 
    mom₂ₜ₋₂ = moments.make_mon_expo(2t-2,M,true)
    @test length(ideal_cons) == length(mom₂ₜ₋₂)*sum(M .== 0.0)

    ξₜⁿⁿⁱᵈ = nn_model.modelξₜⁿⁿⁱᵈ(M,t)
    res = nn_model.computeξₜⁿⁿ(ξₜⁿⁿⁱᵈ)
end

# sparse

@testset "simple sparse t = 1" begin
    M =[1 1 0 0
        1 0 1 0
        0 1 0 1
        0 0 1 1 ] .+ 0.0
    t = 1
    ξₜⁿⁿˢᵖ = nn_model.modelξₜⁿⁿˢᵖ(M,t)
    res = nn_model.computeξₜⁿⁿ(ξₜⁿⁿˢᵖ)
    @test round(objective_value(res), digits=2) == 4.0
end

@testset "simple sparse t = 2" begin
    M = [1 1 0 0
        1 0 1 0
        0 1 0 1
        0 0 1 1 ] .+ 0.0

    t = 2
    ξₜⁿⁿˢᵖ = nn_model.modelξₜⁿⁿˢᵖ(M,t)
    res = nn_model.computeξₜⁿⁿ(ξₜⁿⁿˢᵖ)
    @test round(objective_value(res), digits=2) == 4.0
end

## Constraints
@testset "Ideal constraint sizes" begin
    m,n = rand(1:10), rand(1:10)
    M = nn_matrices.make_NN_mat(m,n)
    nze = moments.get_nonzero_entries(M, false)
    ze = moments.get_zero_entries(M, false)
    model = Model()
    t = 2
    @variable(model, Y[mom.make_mon_expo(m+n,2*t)] ) # Define variables in the moment matrix.
    Lxx = nn_model.get_Lxx(M,t,Y)
    sec_ord_cons = nn_model.get_sec_ord_cons(M,Y)
    @test length(sec_ord_cons) == length(nze)

    loc1_cons   = nn_model.get_loc1_cons(M,t,Y)
    @test length(loc1_cons) == m+n
    for c in loc1_cons
        @test size(c)[1] == length(mom.make_mon_expo(t-1,M,true))
    end
    loc2_cons   = nn_model.get_loc2_cons(M,t,Y)
    @test length(loc2_cons) == length(nze)
    for c in loc2_cons
        @test size(c)[1] == length(moments.make_mon_expo(t-1,M,true))
    end
    ideal_cons  = nn_model.get_ideal_cons(M,t,Y)
    @test length(ideal_cons) == length(ze)*length(mom.make_mon_expo(2t-2, M, true))
    dagger_cons = nn_model.get_dagger_cons(M,t,Y,false)
    @test length(dagger_cons) == length(nze)
    ddagger_cons = nn_model.get_dagger_cons(M,t,Y,true)
    @test length(ddagger_cons) == 1+m+n+length(nze)
end


@testset "sparse constraint sizes" begin
    M = [1 1 0 0
         1 0 1 0
         0 1 0 1
         0 0 1 1 ] .+ 0.0
    t = 1
    m,n = size(M)
    nze = moments.get_nonzero_entries(M, false)
    ze = moments.get_zero_entries(M, false)
    mc = moments.get_maximal_cliques(M, true)
    p = length(mc)
    model = Model()
    spar_inds = nn_model.make_spar_inds(mc,t) 
    @test length(spar_inds) == sum([binomial(length(c)+2t,2t) for c in mc])  # ∑ₖᵖ (|Vₖ + t| choose t) many variables
    @variable(model, Y[spar_inds]) 
    L_kxxV_k = nn_model.get_L_kxxV_k(M,Y,t)
    @test length(L_kxxV_k) == length(mc)
    @test all([size(L_kxxV_k[k])[1] == binomial(length(mc[k])+t,t) for k in 1:p]) # check lengths of vectors
    # @test all([b ⊆ mc[k] for k ∈ 1:p for b ∈ unique(map(a -> findall(a .> 0), Iᵏs[k]))]) # Check if the monomial clques are limited to the resp. cliques

    sparse_sec_ord_cons = nn_model.get_sparse_sec_ord_cons(M,Y)
    @test length(sparse_sec_ord_cons) == length(nze)
    sparse_loc1_cons = nn_model.get_sparse_loc1_cons(M,t,Y)
    @test length(sparse_loc1_cons) == sum(length.(mc))
    sparse_loc2_cons = nn_model.get_sparse_loc2_cons(M,t,Y)
    [M[i,j-m] for (i,j) ∈ moments.get_edges(mc[k],M) for k in 1:p ]  
    
end

@testset "ξₜⁿⁿⁱᵈ < ξₜⁿⁿˢᵖ t = 1" begin
    for i in 1:10
         M = nn_matrices.make_NN_mat(m,n,div(n*m,2))
         t = 1
         ξₜⁿⁿⁱᵈ = nn_model.modelξₜⁿⁿⁱᵈ(M,t) ;
         res_id = nn_model.computeξₜⁿⁿ(ξₜⁿⁿⁱᵈ) ;

         ξₜⁿⁿˢᵖ = nn_model.modelξₜⁿⁿˢᵖ(M,t) ;
         res_sp = nn_model.computeξₜⁿⁿ(ξₜⁿⁿˢᵖ) ;

         @test objective_value(res_id) < objective_value(res_sp)
    end
end

@testset "ξₜⁿⁿⁱᵈ < ξₜⁿⁿˢᵖ t = 2" begin
    for i in 1:10
         M = nn_matrices.make_NN_mat(m,n,div(n*m,2))
         t = 2
         ξₜⁿⁿⁱᵈ = nn_model.modelξₜⁿⁿⁱᵈ(M,t) ;
         res_id = nn_model.computeξₜⁿⁿ(ξₜⁿⁿⁱᵈ) ;

         ξₜⁿⁿˢᵖ = nn_model.modelξₜⁿⁿˢᵖ(M,t) ;
         res_sp = nn_model.computeξₜⁿⁿ(ξₜⁿⁿˢᵖ) ;

         @test objective_value(res_id) < objective_value(res_sp)
    end
end






###



