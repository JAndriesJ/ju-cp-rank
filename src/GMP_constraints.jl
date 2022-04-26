module GMP_constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  make_spar_inds, 
        make_spar_sec_ord_mom_con,
        make_spar_mom_mat_con,
        make_spar_pos_con,
        make_spar_loc_con,
        make_spar_G_con,
        make_obj,
        make_dens_inds, 
        make_dens_sec_ord_mom_con,
        make_dens_mom_mat_con,
        make_dens_pos_con,
        make_dens_loc_con,
        make_dens_G_con,
        make_dens_ideal_con


############################################### Dense GMP constraints#########################################
make_dens_inds(A,t) = mom.make_mon_expo(size(A)[1],2*t)

"""Lᵧ(xxᵀ) = A"""
make_dens_sec_ord_mom_con(A,Y) =  map(α -> Y[α], mom.make_mon_expo(size(A)[1],(1,1),isle=false))

"""Lᵧ([x][x]ᵀ) ⪰ 0""" 
make_dens_mom_mat_con(A,t,Y) = map(α -> Y[α], mom.make_mon_expo(size(A)[1],(t,t)))

"""Lᵧ(xᵢ [x]ₜ₋₁[x]ᵀₜ₋₁ ) for i∈ [n]"""
function make_dens_pos_con(A,t,Y) 
    n = size(A)[1]
    xxₜ₋₁ = mom.make_mon_expo(n,(t-1,t-1))
    return [ map(α -> Y[α+mom.eᵢ(n,i)], xxₜ₋₁)  for i ∈ 1:n]
end

"""Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , Aᵢⱼ ̸= 0, k ∈ [p]"""
function make_dens_loc_con(A,t,Y)
    n = size(A)[1]
    nzes = mom.get_nonzero_entries(A)
    xxₜ₋₁ = mom.make_mon_expo(n,(t-1,t-1))
    lt(nze) = A[nze...]*map(α -> Y[α],xxₜ₋₁)
    rt(nze) = map(α -> Y[α+eᵢ(n,nze[1],nze[2])],xxₜ₋₁)
    return [ lt(nze)-rt(nze)  for nze ∈ nzes]
end

"""Lᵧ([x]ₜ₋₁[x]ₜ₋₁ᵀ ⊗ (A-xxᵀ)) ⪰ 0"""
function make_dens_G_con(A,t,Y)
    n = size(A)[1]
    xxᵀ = mom.make_mon_expo(n,(1,1),isle=false) 
    xxₜ₋₁ = mom.make_mon_expo(n,(t-1,t-1))
    return kron(map(α -> Y[α], xxₜ₋₁), A) - map(α -> Y[α], moments.expo_kron(xxₜ₋₁,xxᵀ))
end
 

"""M(xᵢxⱼL) = 0  for all {i,j} s.t. Mᵢⱼ = 0"""
function make_dens_ideal_con(A,t,Y) 
    n = size(A)[1] 
    momₜ₋₁ = mom.make_mon_expo(t-1,A)
    eᵢⱼs = map( e -> mom.eᵢ(n,e[1],e[2]), mom.get_zero_entries(A))
    return map(eᵢⱼ -> [Y[α+β+eᵢⱼ] for α ∈ momₜ₋₁, β ∈ momₜ₋₁ ], eᵢⱼs)
end


############################################### Sparse GMP constraints#########################################

make_spar_inds(mc,t) = [ [k,m] for  k ∈ 1:length(mc) for m in mom.make_mon_expo(length(mc[k]),2*t) ]

""" ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]"""
function make_spar_sec_ord_mom_con(M,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    Iᵏs₁ = moments.get_monomial_cliques((1,1),M)[1:end-1]
    K = [map(α -> (α ∈ Iᵏs₁[k]) ? Y[[k,α[mc[k]]]] : 0 ,xx₁) for k ∈ 1:p]
    return [(sum(K)[ij...],M[ij...])  for ij ∈ cat([[(i,j) for j ∈ i:n] for i ∈ 1:n]...,dims = 1) if M[ij...] != 0]
end

"""Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]""" 
function make_spar_mom_mat_con(M,t,Y)
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = moments.get_monomial_cliques((t,t),M)[1:end-1]
    Iᵏs_cut = [map(α -> [k,α[mc[k]]], Iᵏs[k]) for k ∈ 1:p]
    return [α_to_Lxᵅ(Y, Iᵏ) for Iᵏ ∈ Iᵏs_cut]
end

"""Lᵧₖ(xᵢ [x]ₜ₋₁ᵥₖ [x]ᵀₜ₋₁ᵥₖ ) for i∈ Vₖ and k ∈ [p]"""
function make_spar_pos_con(M,t,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = get_monomial_cliques(t-1,M)
    arra = [
            [
             [Y[[k, (α + β + eᵢ(n,i))[mc[k]] ]]  # define each monomial (and clique index)
                    for α ∈ Iᵏs[k], β∈ Iᵏs[k]]  # iterate over monomial exponent ass. to a clique
                        for i ∈ mc[k]]    # iterate over each element of a clique
                            for k ∈ 1:p]  # iterate over each clique
    return cat(arra...,dims=1)                        
end

"""Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]"""
function make_spar_loc_con(M,t,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = get_monomial_cliques(t-1,M)
    arra = [
            [
                [temp_fun(M,n,mc,nze,k,Y,α,β) for α ∈ Iᵏs[k], β∈ Iᵏs[k]]  # iterate over monomial exponent ass. to a clique
                        for nze ∈ nz(mc[k])]     # iterate over each pair of element of the clique
                            for k ∈ 1:p]        # iterate over each clique
    return cat(arra...,dims=1) 
end
nz(m) = unique([sort([i,j]) for i in m, j in m])
temp_fun(M,n,mc,nze,k,Y,α,β) = M[nze...]*Y[ [k,(α+β)[mc[k]] ]] - Y[ [k,(α+β + eᵢ(n,nze[1],nze[2]))[mc[k]] ]] 

"""Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]"""
function make_spar_G_con(A,t,Y)
    n = size(A)[1]
    mc = get_maximal_cliques(A)
    p = length(mc)
    Iᵏs = get_monomial_cliques((t-1,t-1),A)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    # lt(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A[mc[k],mc[k]])       # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A|ᵥₖ)
    # rt(k) = map(α -> Y[[k,α[mc[k]]]], expo_kron(Iᵏs[k],xx₁[mc[k],mc[k]])) # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ)
    lt(k) = kron(map(α -> Y[[k,α[mc[k]]]], Iᵏs[k]), A)                                  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ A)
    rt(k) = map(α ->  isinmc(α,mc[k]) ? Y[[k,α[mc[k]]]] : 0, expo_kron(Iᵏs[k], xx₁))    # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|ᵥₖ₌₀  )
    return [lt(k) - rt(k) for k in 1:p]
    # [get_t(A,xx₁,mc,k,Iᵏs,Y) for k in 1:p] 
end
isinmc(α,mck) = all([j ∈ mck for j ∈ findall(α .> 0)])
# get_t(M,xx₁,mc,k,Iᵏs,Y) = get_lt(M,mc,k,Iᵏs,Y) - get_rt(xx₁,mc,k,Iᵏs,Y)                # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ)
# get_rt(xx₁,mc,k,Iᵏs,Y) =  α_to_Lxᵅ(Y,moments.expo_kron(Iᵏs[k],xx₁[mc[k],mc[k]]),mc[k])  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ xxᵀ|_vₖ)
# get_lt(M,mc,k,Iᵏs,Y) = kron(α_to_Lxᵅ(Y, Iᵏs[k],mc[k]), M[mc[k],mc[k]])                  # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⊗ A|_vₖ               

"""∑ Lᵧₖ(1)"""
make_obj(mc,Y) = sum([Y[[k,zeros(Int,length(mc[k]))]] for k in 1:length(mc)])

"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Y, index_array, mck) =  map(γ -> Y[[k,γ[mck]]], index_array)
α_to_Lxᵅ(Lx, index_array) = map(α -> Lx[α], index_array)

end