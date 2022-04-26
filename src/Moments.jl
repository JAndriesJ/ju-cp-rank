module moments
using Test
using Graphs
using LinearAlgebra

export  eᵢ,
        make_mon_expo,
        get_maximal_cliques,
        get_monomial_cliques,
        expo_kron,
        get_mat_support,
        get_mom_cliq_supp,
        get_mom_mat_supp,
        get_ten_con_supp,
        run_tests
       
        

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]
function eᵢ(n::Int,i::Int,j::Int) 
    if i == j
        return [k == i ? 2 : 0 for k in 1:n]
    else
        return [k ∈ [i,j] ? 1 : 0 for k in 1:n]
    end
end

"""[x]≦ₜ := [xᵅ for all α ∈ ℕⁿₜ] or [x]₌ₜ  := [xᵅ for all α ∈ ℕⁿ≤ₜ]"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    @assert typeof(t) == Int64
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp,M_vec) : M_vec)
end
"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true) = make_mon_expo(make_mon_expo(n,t[1]; isle=isle), make_mon_expo(n,t[2]; isle=isle))
make_mon_expo(mom::Vector{Vector{Int64}}) = make_mon_expo(mom,mom) 
make_mon_expo(mom₁::Vector{Vector{Int64}},mom₂::Vector{Vector{Int64}}) = [α+β for α ∈ mom₁, β ∈ mom₂]
make_mon_expo(t::Int,M::Matrix{Float64}) = unique(cat(get_monomial_cliques(t,M)[1:end-1]...,dims=1))
make_mon_expo(t::Tuple{Int,Int},M) = make_mon_expo(make_mon_expo(t[1],M),make_mon_expo(t[2],M))
make_mon_expo(t::Int,M::Matrix{Float64},j) = unique(cat(get_monomial_cliques(t,M,j)[1:end-1]...,dims=1))
make_mon_expo(t::Tuple{Int,Int},M,j) = make_mon_expo(make_mon_expo(t[1],M,j),make_mon_expo(t[2],M,j))

## Sparsity 
### Cliques
get_maximal_cliques(M) = Graphs.maximal_cliques(Graph(M .> 0))

"""[x]≦ₜ ᵥₖ[x]ᵀ≦ₜ ᵥₖ """
function get_monomial_cliques(n::Int,t::Int,mc::Vector{Vector{Int64}}) 
    Iks = [get_mon_clique(n,t,c) for c in mc]
    Iks_comp = setdiff(make_mon_expo(n,t), union(Iks...))
    return [Iks..., Iks_comp]
end
get_monomial_cliques(t::Int,M::Matrix{Float64}) = get_monomial_cliques(size(M)[1],t,get_maximal_cliques(M))
get_monomial_cliques(t::Tuple{Int,Int},M::Matrix{Float64}) = [make_mon_expo(m,m) for m in get_monomial_cliques(t[1],M)]
function get_monomial_cliques(t::Int,M::Matrix{Float64},j::Int)
    n = size(M)[1]
    @assert j ∈ [1:n...]
    Ik_s = get_monomial_cliques(t,M)
    supp_of_Iks = get_supp_of_Iks(Ik_s)
    Ik_s_j  = push!(j .∈ supp_of_Iks, 0)
    Ik_s_nj = push!(j .∉ supp_of_Iks, 1)
    return cat(Ik_s[Ik_s_j],[union(Ik_s[Ik_s_nj]...)],dims=1)
end
get_monomial_cliques(t::Tuple{Int,Int},M::Matrix{Float64},j::Int) = [make_mon_expo(m,m) for m in get_monomial_cliques(t[1],M,j)]

get_supp_of_Iks(Ik_s) = [union(map(i->findall(i .> 0),I)...) for I in Ik_s[1:end-1] ]
get_mon_clique(n,t,c) = map(v->embed(v,c,n), make_mon_expo(length(c),t))
embed(v,α,n) = [i ∈ α ? popfirst!(v) : 0  for i in 1:n]


### supports
get_supp_mat(M) = (M .> 0.0) .+ 0
get_supp_mom_mat(t,M) = (sum(get_supp_mom_cliq(t,M)) .> 0) .+ 0
get_supp_ten_con(t,M) = kron(get_supp_mat(M), get_mom_mat_supp(t,M)) 
function get_supp_mom_cliq(t,M)
    MCs = get_monomial_cliques(t,M)
    Yᵏs = [[[α,β] for α ∈ m, β ∈ m] for m in MCs]
    mvₜ = make_mon_expo(t,M)
    return [[ any([α,β] ∈ Yᵏ) ? 1 : 0 for α ∈ mvₜ , β ∈ mvₜ] for Yᵏ in Yᵏs ]
end

## Utilities
"""A ∈ (ℕⁿ)ᵃˣᵇ, B ∈ (ℕⁿ)ᶜˣᵈ --> D ∈ (ℕⁿ)ᵃᶜˣᵇᵈ : D₍ᵢⱼ,ₖₕ₎ = Aᵢₖ + Bⱼₕ"""
function expo_kron(A,B)
    n₁,n₂ = size(A)
    D = [B + repeat( [A[i,j]] , inner = (1,1), outer = size(B)) for i in 1:n₁ , j in 1:n₂ ]
    return cat([cat(D[i,:]...,dims=2) for i in 1:n₁]...,dims=1)
end

get_zero_entries(M) = [ (i,j) for i in 1:(size(M)[1]-1) for j in (i+1):size(M)[1] if M[i,j] == 0]
get_nonzero_entries(M) = [(i,j) for i in 1:size(M)[1] for j in i:size(M)[1] if M[i,j] != 0]



### Tests
function run_tests()
    @testset "standerd basis vector" begin
        for n in 6:13
            for k in 1:n
                e = eᵢ(n,k)
                @test length(e) == n
                @test maximum(e) == 1
                @test minimum(e) == 0
                @test sum(e) == 1
            end
        end

        # I = rand(1:9,4)
        # e = eᵢ(maximum(I),I)
        # @test length(e) == maximum(I)
        # @test sum(e .== 0) == (maximum(I) - 4)
    end

    @testset "make_mon_expo" begin
        for n ∈ 3:9, t ∈ 0:4
            @test size.(make_mon_expo(n,t,isle=true))[1] == (n,)
            @test length(make_mon_expo(n,t,isle=true)) == binomial(n+t,t)
        end

        MonBase = make_mon_expo(2,3,isle=false)
        @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]

        MonBase = make_mon_expo(2,4,isle=false)
        @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]

        MB = make_mon_expo(2,(3,3),isle=true)
        @test size(MB) == (binomial(2+3,3),binomial(2+3,3))

        @test make_mon_expo(make_mon_expo(2,3)) == make_mon_expo(2,(3,3))
    end

    @testset "expo_kron" begin
        n = 4
        t = (2,3)
        A = make_mon_expo(n,t)
        B = make_mon_expo(n,t)
        AOXB = expo_kron(A,B)
        @test size(AOXB) == size(A) .* size(B) 
    end

    @testset "Misc" begin
        M = [1 0 ; 0 1]
        println(get_zero_entries(M))
        println(get_nonzero_entries(M))
        @test get_zero_entries(M) == [(1, 2)]
        @test get_nonzero_entries(M) == [(1, 1), (2, 2)]

        @test get_supp_mat(M) == M
        # get_supp_mom_mat(t,M) 
        # get_supp_ten_con(t,M) 

    end

end

end

