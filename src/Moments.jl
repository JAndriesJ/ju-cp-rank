module moments
using Test
using Graphs

export  eᵢ,
        eᵢⱼ,
        make_mon_expo,
        get_monomial_cliques,
        run_tests

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]
eᵢⱼ(n::Int,i::Int,j::Int) = [k ∈ [i,j] ? 1 : 0 for k in 1:n]

"""[x]≦ₜ := [xᵅ for all α ∈ ℕⁿₜ] or [x]₌ₜ  := [xᵅ for all α ∈ ℕⁿ≤ₜ]"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    @assert typeof(t) == Int64
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp,M_vec) : M_vec)
end

"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
function make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true)
    M_vec1      = make_mon_expo(n,t[1]; isle=isle)
    M_vec2      = make_mon_expo(n,t[2]; isle=isle)
    return [mi+mj for mi in M_vec1, mj in M_vec2]
end

make_mon_expo(mom::Vector{Vector{Int64}}) = make_mon_expo(mom,mom) 
make_mon_expo(mom₁::Vector{Vector{Int64}},mom₂::Vector{Vector{Int64}}) = [α+β for α ∈ mom₁, β ∈ mom₂]

### Sparsity 

get_zero_entries(M) = [ (i,j) for i in 1:(size(M)[1]-1) for j in (i+1):size(M)[1] if M[i,j] == 0]
get_nonzero_entries(M) = [(i,j) for i in 1:(size(M)[1]-1) for j in i:size(M)[1] if M[i,j] != 0]
get_maximal_cliques(M) = Graphs.maximal_cliques(Graph(M .> 0))
#---
# function get_monomial_cliques(n,t,M,j) 
#     @assert j ∈ [1:n...]
#     Iks = get_monomial_cliques(n,t,M)
#     [ I in Iks]

# end


function get_monomial_cliques(n,t,mc::Vector{Vector{Int64}}) 
    Iks = [get_mon_clique(n,t,c) for c in mc]
    Iks[1] = unique([intersect(Iks...)..., Iks[1]...]) # put the intersetion first
    mon = make_mon_expo(n,t)
    Iks_comp = setdiff(mon, union(Iks...))
    return [Iks..., Iks_comp]
end
get_monomial_cliques(n,t,M::Matrix{Float64}) = get_monomial_cliques(n,t,get_maximal_cliques(M))
function get_monomial_cliques(n,t::Int,M::Matrix{Float64},j::Int)
    @assert j ∈ [1:n...]
    Ik_s = get_monomial_cliques(n,t,M)
    supp_of_Iks = get_supp_of_Iks(Ik_s)
    Ik_s_j  = push!(j .∈ supp_of_Iks, 0)
    Ik_s_nj = push!(j .∉ supp_of_Iks, 1)
    return cat(Ik_s[Ik_s_j],[union(Ik_s[Ik_s_nj]...)],dims=1)
end
get_supp_of_Iks(Ik_s) = [union(map(i->findall(i .> 0),I)...) for I in Ik_s[1:end-1] ]
get_mon_clique(n,t,c) = map(v->embed(v,c,n), make_mon_expo(length(c),t))
embed(v,α,n) = [i ∈ α ? popfirst!(v) : 0  for i in 1:n]
#---
make_mon_expo(n,t,mc::Vector{Vector{Int64}}) = unique(cat(get_monomial_cliques(n,t,mc)[1:end-1]...,dims=1))
make_mon_expo(n,t::Int,M::Matrix{Float64}) = unique(cat(get_monomial_cliques(n,t,M)[1:end-1]...,dims=1))
make_mon_expo(n,t::Tuple{Int,Int},M) = make_mon_expo(make_mon_expo(n,t[1],M),make_mon_expo(n,t[2],M))
#---
make_mon_expo(n,t::Int,M::Matrix{Float64},j) = unique(cat(get_monomial_cliques(n,t,M,j)[1:end-1]...,dims=1))
make_mon_expo(n,t::Tuple{Int,Int},M,j) = make_mon_expo(make_mon_expo(n,t[1],M),make_mon_expo(n,t[2],M))


function order_monomials(mon, mon_cliques)
    k = length(mon_cliques) - 1
    mon_cli_index  = map(m -> findall([m in mc for mc in mon_cliques]), mon)
    sp = sortperm(mon_cli_index)
    intersection = findall(map(c -> c == 1:k, mon_cli_index))
    sp = setdiff(sp,intersection)
    return cat(intersection,sp,dims=1)
end


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
end

end
