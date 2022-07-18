module constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  make_loc_con, # Dense
        make_G_con,
        make_ideal_con, 
        make_dag_con, 
        make_xx_con,
        make_spar_inds, # Sparse
        make_spar_sec_ord_mom_con, 
        make_spar_mom_mat_con, 
        make_spar_dag_con, 
        make_spar_xx_con,  
        make_spar_loc_con, 
        make_spar_G_con, 
        make_obj, 
        α_to_Lxᵅ
 
        


############################################# Sparse constraints#########################################

make_spar_inds(mc,t) = [ [k,m] for  k ∈ 1:length(mc) for m in mom.make_mon_expo(length(mc[k]), 2*t) ]

""" ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]"""
function make_spar_sec_ord_mom_con(M,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    Iᵏs₁ = mom.get_monomial_cliques((1,1), M)[1:end-1]
    K = sum([map(α -> (α ∈ Iᵏs₁[k]) ? Y[[k,α[mc[k]]]] : 0, xx₁) for k ∈ 1:p])
    # ut = [(i,j) for i ∈ 1:n for j ∈ i:n if M[i,j] != 0 ] #cat([[(i,j) for j ∈ i:n] for i ∈ 1:n]...,dims = 1)
    nze = mom.get_nonzero_entries(M)
    return [(K[ij...], M[ij...])  for ij ∈ nze]
end

"""Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]""" 
function make_spar_mom_mat_con(M,t,Y)
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = mom.get_monomial_cliques((t,t),M)[1:end-1]
    Iᵏs_cut = [map(a -> [k, a[mc[k]] ], Iᵏs[k]) for k ∈ 1:p]
    return [α_to_Lxᵅ(Y, Iᵏ) for Iᵏ ∈ Iᵏs_cut]
end

"""
L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]
L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
"""
function make_spar_dag_con(A,t,Y, dd=false)
    n = size(A)[1]
    mc = moments.get_maximal_cliques(A)
    p = length(mc)
    moms₂ₜ₋₂ = moments.get_monomial_cliques(2*t-2,A)[1:end-1]
    moms₂ₜ   = moments.get_monomial_cliques(2*t,A)[1:end-1]
    # L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
    oddone = [[Y[[k, a[mc[k]]]] for a in moms₂ₜ[k]] for k in 1:p]
    # L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]   
    ondiag = [[sqrt(A[i,i])*Y[[k, (mom.eᵢ(n,i)+a)[mc[k]]]] - Y[[k, (2*mom.eᵢ(n,i) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for i ∈ mc[k]]
    # L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
    offdiag = [[A[i,j]*Y[[k, a[mc[k]]]] - Y[[k, (mom.eᵢ(n,i) + mom.eᵢ(n,j) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for i ∈ mc[k] for j ∈ setdiff(Set(mc[k]),Set([i]))]
    # warning("This could be wrong")
    return  dd ? [oddone...,offdiag...,  ondiag...] : offdiag 
end

"""
L(xᵢxⱼu) for u ∈ [x]ᵥₖ_₂ₜ₋₂ ; i,j ∈ mc[k], k ∈ [p] 
"""
function make_spar_xx_con(A,t,Y)
    n = size(A)[1]
    mc = moments.get_maximal_cliques(A)
    p = length(mc)
    momsₜ₋₁ = moments.get_monomial_cliques(t-1,A)
    # L(xᵢxⱼu) for u ∈ [x]ᵥₖ_₂ₜ₋₂ ; i,j ∈ mc[k],  k ∈ [p]  
    return [
            [ Y[[k,(a + b + mom.eᵢ(n,i) + mom.eᵢ(n,j))[mc[k]]]]
               for a in momsₜ₋₁[k], b ∈ momsₜ₋₁[k] 
                                                 ] for k ∈ 1:p for i ∈ mc[k] for j ∈ setdiff(Set(mc[k]),Set([i]))  ]
end

""" Lₖ((√Aᵢᵢxᵢ-xᵢ²)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ)  ⪰ 0, i ∈ Vₖ
    Lₖ((Aᵢⱼ-xᵢxⱼ)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0, i,j ∈ Vₖ, k ∈ [p] """
function make_spar_loc_con(M,t,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = get_monomial_cliques(t-1,M)

    term(e,k,a,b)     = e[1]==e[2] ? d_term(k,e,a,b) : od_term(k,e,a,b)
    d_term(k,e,a,b)   = sqrt(M[e...])*Y[ [k, (a+b+eᵢ(n,e[1]) )[mc[k]]] ] - Y[ [k, (a+b + eᵢ(n,e[1],e[2]))[mc[k]]]] # Lₖ((√Mᵢᵢxᵢ-xᵢ²)xᵅ)
    od_term(k,e,a,β)  =      M[e...] *Y[ [k, (a+β            )[mc[k]]] ] - Y[ [k, (a+β + eᵢ(n,e[1],e[2]))[mc[k]]]] # Lₖ((Aᵢⱼ-xᵢxⱼ)xᵅ)

    arr = [[[term(e,k,a,b)  for a ∈ Iᵏs[k], b ∈ Iᵏs[k]]      # iterate over monomial exponent ass. to a clique
                            for e ∈ get_edges(sort(mc[k]))]  # iterate over each pair of element of the clique
                            for k ∈ 1:p]                     # iterate over each clique
    return cat(arr...,dims=1) 
end
get_edges(c) = [[c[i],c[j]] for i ∈ 1:length(c) for j ∈ i:length(c)]


"""Lᵧₖ( [x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]"""
function make_spar_G_con(A,t,Y,isWeak=false)
    n = size(A)[1]
    mc = get_maximal_cliques(A)
    p = length(mc)
    Iᵏs = mom.get_monomial_cliques((t-1,t-1),A)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    
    if isWeak
        lt_w(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A[mc[k],mc[k]])                       # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A|ᵥₖ)
        rt_w(k) = map(α -> Y[[k,α[mc[k]]]], mom.expo_kron(Iᵏs[k], xx₁[mc[k],mc[k]]))            # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ)
        return [lt_w(k) - rt_w(k) for k in 1:p]
    else
        lt(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A)                                       # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A)
        rt(k) = map(α ->  isinmc(α, mc[k]) ? Y[[k,α[mc[k]]]] : 0, mom.expo_kron(Iᵏs[k], xx₁))    # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ₌₀  )
        return [lt(k) - rt(k) for k in 1:p]
    end
end
isinmc(α,mck) = all([j ∈ mck for j ∈ findall(α .> 0)])

"""∑ₖ Lₖ(1)"""
make_obj(mc,Y) = sum([Y[[k,zeros(Int,length(mc[k]))]] for k in 1:length(mc)])

############################################### Utility #########################################
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Y, index_array) = map(α -> Y[α], index_array)
α_to_Lxᵅ(Y, index_array, mck) =  map(γ -> Y[[k,γ[mck]]], index_array)

end