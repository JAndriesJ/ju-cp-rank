module moments
using Test
using Graphs
 
export  eᵢ,
        make_mon_expo,
        get_maximal_cliques,
        get_monomial_cliques,
        expo_kron

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]
eᵢ(n::Int,i::Int,j::Int) = i == j ? [k == i ? 2 : 0 for k in 1:n] : [k ∈ [i,j] ? 1 : 0 for k in 1:n] 
#make_mon_expo-------------------------------------------------------------------------------------------------------------------------------------
"""[x]≦ₜ := [xᵅ for all α ∈ ℕⁿₜ] or [x]₌ₜ  := [xᵅ for all α ∈ ℕⁿ≤ₜ]"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp, M_vec) : M_vec)
end
"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true) = make_mon_expo(make_mon_expo(n,t[1]; isle=isle), make_mon_expo(n,t[2]; isle=isle))
make_mon_expo(mom₁::Vector{Vector{Int64}},mom₂::Vector{Vector{Int64}}) = [a+b for a ∈ mom₁, b ∈ mom₂]
make_mon_expo(mom::Vector{Vector{Int64}}) = make_mon_expo(mom, mom) 
# sparse--------------------------------------------------------------------------------------------------------------------------------------------
"""[x(Vₖ)]ₜ[x(Vₖ)]ᵀₜ""" 
# make_mon_expo(t::Int,M::Matrix{Float64}) = unique(cat(get_monomial_cliques(t,M)[1:end-1]...,dims=1))
make_mon_expo(t::Int,M::Matrix{Float64},j::Int) = t != 0 ? unique(cat(get_monomial_cliques(t,M,j)[1:end-1]...,dims=1)) : get_monomial_cliques(t,M,j)[1]
function make_mon_expo(t, M, is_nn=false; isle::Bool = true)
    m,n = size(M)
    if is_nn
        mom₂ₜ = make_mon_expo(m+n,t, isle=isle)
        mc = moments.get_maximal_cliques(M, is_nn)
        push!(mc,[1:m...],[m+1:m+n...])
        return [m for m in mom₂ₜ if check_sup(m,mc)]
    else
        mom₂ₜ = make_mon_expo(n,t, isle=isle)
        return mom₂ₜ
    end
end
check_sup(mom,mc) = any(map(c -> findall(mom .> 0) ⊆ c, mc))
make_mon_expo(t::Int,M::Matrix{Float64},j::Int,isnn=false) = t != 0 ? unique(cat(get_monomial_cliques(t,M,j,isnn)[1:end-1]...,dims=1)) : get_monomial_cliques(t,M,j,isnn)[1]
#get_monomial_cliques-------------------------------------------------------------------------------------------------------------------------------
"""[x(Vₖ)]ₜ[x(Vₖ)]ᵀₜ """
function get_monomial_cliques(n::Int,t::Int,mc::Vector{Vector{Int64}}) 
    Iks = [get_mon_clique(n,t,c) for c in mc]
    Iks_comp = setdiff(make_mon_expo(n,t), union(Iks...))
    return [Iks..., Iks_comp]
end
get_monomial_cliques(t::Tuple{Int,Int},M::Matrix{Float64},isnn=false::Bool) = [make_mon_expo(m,m) for m in get_monomial_cliques(t[1],M,isnn)]
function get_monomial_cliques(t::Int,M::Matrix{Float64},isnn=false)
    if isnn
        return get_monomial_cliques(sum(size(M)),t, get_maximal_cliques(M,isnn))
    else
        return get_monomial_cliques(size(M)[1],t, get_maximal_cliques(M))
    end
end
function get_monomial_cliques(t::Int,M::Matrix{Float64},j::Int,isnn=false)
    m,n = size(M)
    if t==0 
         return  isnn ? [[zeros(Int,m+n)]]  : [[zeros(Int,n)]]
    else
        Ik_s = get_monomial_cliques(t,M,isnn) 
        supp_of_Iks = get_supp_of_Iks(Ik_s)
        Ik_s_j  = push!(j .∈ supp_of_Iks, 0)
        Ik_s_nj = push!(j .∉ supp_of_Iks, 1)
        return cat(Ik_s[Ik_s_j],[union(Ik_s[Ik_s_nj]...)],dims=1)
    end
end
get_mon_clique(n,t,c) = map(v->embed(v,c,n), make_mon_expo(length(c),t))
embed(v,α,n) = [i ∈ α ? popfirst!(v) : 0  for i in 1:n]
get_supp_of_Iks(Ik_s) = [union(map(i->findall(i .> 0),I)...) for I in Ik_s[1:end-1] ]
#cliques--------------------------------------------------------------------------------------------------------------------------------------------
function get_maximal_cliques(M, isnn=false)
    if isnn
         m,n = size(M)
         Gₘ_adj_mat = make_Gₘ_adj_mat(M) 
         mc = sort.(Graphs.maximal_cliques(Graph(Gₘ_adj_mat .> 0))) 
         return [c for c ∈ mc if !(c ⊆ 1:m) && !(c ⊆ m+1:m+n)]
    else
        return sort.(Graphs.maximal_cliques(Graph(M .!= 0)))
    end
end
function get_edges(c::Vector{Int64},loops=true::Bool) 
    k = length(c)
    if loops 
        return [[c[i],c[j]] for i ∈ 1:k for j ∈ i:k] 
    else
        return [[c[i],c[j]] for i ∈ 1:k for j ∈ i:k if i != j]
    end
end
function get_edges(c::Vector{Int64},M::Matrix{Float64}) 
    m,n = size(M)
    return [[a,b] for a ∈ c∩[1:m...] for b ∈ c∩[m+1:m+n...]] 
end

function make_Gₘ_adj_mat(M)
    r,c = size(M)
    M_supp = M .> 0.0
    return vcat(hcat(ones(r,r),  M_supp), 
                hcat( M_supp',ones(c,c)))
end
## Utilities-----------------------------------------------------------------------------------------------------------------------------------------
"""A ∈ (ℕⁿ)ᵃˣᵇ, B ∈ (ℕⁿ)ᶜˣᵈ --> D ∈ (ℕⁿ)ᵃᶜˣᵇᵈ : D₍ᵢⱼ,ₖₕ₎ = Aᵢₖ + Bⱼₕ"""
function expo_kron(A,B)
    n₁,n₂ = size(A)
    D = [B + repeat( [A[i,j]] , inner = (1,1), outer = size(B)) for i in 1:n₁ , j in 1:n₂ ]
    return cat([cat(D[i,:]...,dims=2) for i in 1:n₁]...,dims=1)
end
"""get all (i,j) s.t. Mᵢⱼ ≠ 0"""
function get_nonzero_entries(M, is_cp=true)
    m,n = size(M)
    if is_cp
        return [(i,j) for i in 1:m for j in i:n if M[i,j] != 0] 
    else
        return [(i,j) for i in 1:m for j in 1:n if M[i,j] != 0]
    end
end
"""get all (i,j) s.t. Mᵢⱼ = 0"""
get_zero_entries(M, is_cp=true) = get_nonzero_entries(M .== 0, is_cp)

end
