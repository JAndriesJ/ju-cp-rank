module GMP_constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  make_spar_inds, 
        make_spar_sec_ord_mom_con,
        make_spar_mom_mat_con,
        make_spar_loc_con,
        make_spar_G_con,
        make_obj

############################################# Sparse constraints#########################################

make_spar_inds(mc,t) = [ [k,m] for  k ∈ 1:length(mc) for m in mom.make_mon_expo(length(mc[k]),2*t) ]

""" ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]"""
function make_spar_sec_ord_mom_con(M,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    Iᵏs₁ = moments.get_monomial_cliques((1,1), M)[1:end-1]
    K = [map(α -> (α ∈ Iᵏs₁[k]) ? Y[[k,α[mc[k]]]] : 0, xx₁) for k ∈ 1:p]
    ut = cat([[(i,j) for j ∈ i:n] for i ∈ 1:n]...,dims = 1)
    return [(sum(K)[ij...],M[ij...])  for ij ∈ ut if M[ij...] != 0]
end

"""Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]""" 
function make_spar_mom_mat_con(M,t,Y)
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = moments.get_monomial_cliques((t,t),M)[1:end-1]
    Iᵏs_cut = [map(α -> [k, α[mc[k]] ], Iᵏs[k]) for k ∈ 1:p]
    return [α_to_Lxᵅ(Y, Iᵏ) for Iᵏ ∈ Iᵏs_cut]
end


""" Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ Vₖ
    Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p] """
function make_spar_loc_con(M,t,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = get_monomial_cliques(t-1,M)
    arr = [[[term(M,n,mc,e,k,Y,α,β) for α ∈ Iᵏs[k], β∈ Iᵏs[k]]       # iterate over monomial exponent ass. to a clique
                                     for e ∈ clique_ele_pairs(mc[k])] # iterate over each pair of element of the clique
                                     for k ∈ 1:p]                     # iterate over each clique
    return cat(arr...,dims=1) 
end
clique_ele_pairs(m) = unique([sort([i,j]) for i ∈ m, j ∈ m])
term(M,n,mc,e,k,Y,α,β) = e[1] == e[2] ? d_term(M,Y,n,mc,k,e,α,β) : od_term(M,Y,n,mc,k,e,α,β)
d_term(M,Y,n,mc,k,e,α,β)   = sqrt(M[e...])*Y[ [k, ( α+β+eᵢ(n,e[1]) )[mc[k]]] ] - Y[ [k, ( α+β + eᵢ(n,e[1],e[2]) )[mc[k]]] ] # Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)xᵅ)
od_term(M,Y,n,mc,k,e,α,β)  =      M[e...] *Y[ [k, ( α+β            )[mc[k]]] ] - Y[ [k, ( α+β + eᵢ(n,e[1],e[2]))[mc[k]] ]]  # Lᵧₖ((Aᵢⱼ-xᵢxⱼ)xᵅ)


"""Lᵧₖ( [x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]"""
function make_spar_G_con(A,t,Y,isWeak=false)
    n = size(A)[1]
    mc = get_maximal_cliques(A)
    p = length(mc)
    Iᵏs = get_monomial_cliques((t-1,t-1),A)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    
    if isWeak
        lt_w(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A[mc[k],mc[k]])                   # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A|ᵥₖ)
        rt_w(k) = map(α -> Y[[k,α[mc[k]]]], expo_kron(Iᵏs[k], xx₁[mc[k],mc[k]]))            # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ)
        return [lt_w(k) - rt_w(k) for k in 1:p]
    else
        lt(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A)                                  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A)
        rt(k) = map(α ->  isinmc(α, mc[k]) ? Y[[k,α[mc[k]]]] : 0, expo_kron(Iᵏs[k], xx₁))    # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ₌₀  )
        return [lt(k) - rt(k) for k in 1:p]
    end
end
isinmc(α,mck) = all([j ∈ mck for j ∈ findall(α .> 0)])

"""∑ Lᵧₖ(1)"""
make_obj(mc,Y) = sum([Y[[k,zeros(Int,length(mc[k]))]] for k in 1:length(mc)])

"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Y, index_array, mck) =  map(γ -> Y[[k,γ[mck]]], index_array)
α_to_Lxᵅ(Lx, index_array) = map(α -> Lx[α], index_array)

end





# """Lᵧ(xᵢ [x]ₜ₋₁[x]ᵀₜ₋₁ ) for i∈ [n]"""
# function make_dens_pos_con(A,t,Y) 
#     n = size(A)[1]
#     xxₜ₋₁ = mom.make_mon_expo(n,(t-1,t-1))
#     return [ map(α -> Y[α+mom.eᵢ(n,i)], xxₜ₋₁)  for i ∈ 1:n]
# end

# """Lᵧₖ(xᵢ [x]ₜ₋₁ᵥₖ [x]ᵀₜ₋₁ᵥₖ ) for i∈ Vₖ and k ∈ [p]"""
# function make_spar_pos_con(M,t,Y)
#     n = size(M)[1]
#     mc = get_maximal_cliques(M)
#     p = length(mc)
#     Iᵏs = get_monomial_cliques(t-1,M)
#     arra = [
#             [
#              [Y[[k, (α + β + eᵢ(n,i))[mc[k]] ]]  # define each monomial (and clique index)
#                     for α ∈ Iᵏs[k], β∈ Iᵏs[k]]  # iterate over monomial exponent ass. to a clique
#                         for i ∈ mc[k]]    # iterate over each element of a clique
#                             for k ∈ 1:p]  # iterate over each clique
#     return cat(arra...,dims=1)                        
# end

# get_t(M,xx₁,mc,k,Iᵏs,Y) = get_lt(M,mc,k,Iᵏs,Y) - get_rt(xx₁,mc,k,Iᵏs,Y)                # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ)
# get_rt(xx₁,mc,k,Iᵏs,Y) =  α_to_Lxᵅ(Y,moments.expo_kron(Iᵏs[k],xx₁[mc[k],mc[k]]),mc[k])  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|_vₖ)
# get_lt(M,mc,k,Iᵏs,Y) = kron(α_to_Lxᵅ(Y, Iᵏs[k],mc[k]), M[mc[k],mc[k]])                  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⊗ A|_vₖ               
