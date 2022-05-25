module constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  α_to_Lxᵅ,
        make_loc_con, # Dense
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
        make_obj 
 
        
"""
input: A(Data matrix),t(integer),Y(JuMP variable)
output: L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)
"""
function make_dag_con(A,t,Y)
    n = size(A)[1]
    mom₂ₜ   = moments.make_mon_expo(2*t,A)
    mom₂ₜ₋₂ = moments.make_mon_expo(2*t-2,A)
    # L(u) for u ∈ [x]₂ₜ
    oddone = [[ 1*Y[α ] for α in mom₂ₜ]]
    # L((√Aₖₖ xₖ - xₖ²)⋅u) for u ∈ [x]₂ₜ₋₂
    ondiag = [[ sqrt(A[i,i])*Y[a] - Y[2*mom.eᵢ(n,i) + a]  for a in mom₂ₜ₋₂] for i ∈ 1:n ]
    # L((Aₖₕ  - xₖxₕ)⋅u) for u ∈ [x]₂ₜ₋₂
    offdiag = [[ A[i,j]*Y[a] - Y[mom.eᵢ(n,i) + mom.eᵢ(n,j) + a]  for a in mom₂ₜ₋₂] for i ∈ 1:n for j ∈ i+1:n]
    return [oddone..., offdiag..., ondiag...]
end
   
"""
input: A(Data matrix), t(integer), Y(JuMP variable)
output: L(xᵢxⱼu) ⪰ 0 for u ∈ [x]₂ₜ₋₂ and i,j ∈ [n]
"""
function make_xx_con(A,t,Y)
    n = size(A)[1]
    mom₁ₜ₋₁ = mom.make_mon_expo(t-1,A)
    # L(xᵢxⱼu) for g ∈ {1} ∪ Sᶜᵖ_A and i,j ∈ [n]
    return [[ Y[a + b + mom.eᵢ(n,i) + mom.eᵢ(n,j)]  for a in mom₁ₜ₋₁, b ∈ mom₁ₜ₋₁] for i ∈ 1:n for j ∈ i:n if A[i,j] != 0.0]
end

"""
input: A(data array), LMB(moment exponent vector), Y(JuMP variable)
output: dictionary: Keys: (i,j) ∈ [n]²
                    vals:
comment: L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ)
(M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function make_loc_con(A,t,Y)
    n = size(A)[1]
    nze = mom.get_nonzero_entries(A)
    momₜ₋₁  = mom.make_mon_expo(t-1, A)
    # if_diag(k) = [sqrt(A[k,k])*Y[α+β+mom.eᵢ(n,k)] - Y[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    # if_off_diag(k,h) = [A[k,h]*Y[α+β] - Y[α+β+mom.eᵢⱼ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    momₜ₋₁g(j)  = mom.make_mon_expo(t-1,A,j)
    if_diag(k) = [sqrt(A[k,k])*Y[α+β+mom.eᵢ(n,k)] - Y[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁g(k), β ∈ momₜ₋₁g(k)]
    if_off_diag(k,h) = [A[k,h]*Y[α+β] - Y[α+β+mom.eᵢ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    return [e[1]==e[2] ? if_diag(e[1]) : if_off_diag(e[1],e[2]) for e in nze]
end

"""M(G ⊗ L) ⪰ 0 constraints
input: A(data matrix),t(Integer),Y(JuMP variable)
Assumption: G = A-[x]₌₁[x]₌₁ᵀ
output: A⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
"""
function make_G_con(A,t,Y)
    n = size(A)[1]
    LMBexp₁    = mom.make_mon_expo(n,(1,1),isle=false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexpₜ₋₁   = mom.make_mon_expo((t-1,t-1),A)     #exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁      = α_to_Lxᵅ(Y,LMBexpₜ₋₁)    #L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp₁ₜ₋₁  = mom.expo_kron(LMBexp₁,LMBexpₜ₋₁)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB₁ₜ₋₁     = α_to_Lxᵅ(Y,LMBexp₁ₜ₋₁)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB₁ₜ₋₁             # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    return G_con
end

"""L(xᵢxⱼ[x]₂ₜ₋₂) = 0  for all {i,j} s.t. Mᵢⱼ = 0"""
function make_ideal_con(A,t,Y) 
    n = size(A)[1] 
    mom₂ₜ₋₂ = mom.make_mon_expo(2t-2, A)
    eᵢⱼs = map( e -> mom.eᵢ(n,e[1],e[2]), mom.get_zero_entries(A))
    return map(eᵢⱼ -> [Y[ α+eᵢⱼ ] for α ∈ mom₂ₜ₋₂ ], eᵢⱼs)
end


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

"""
L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]
L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
"""
function make_spar_dag_con(A,t,Y)
    n = size(A)[1]
    mc = moments.get_maximal_cliques(A)
    p = length(mc)
    moms₂ₜ₋₂ = moments.get_monomial_cliques(2*t-2,A)[1:end-1]
    moms₂ₜ   = moments.get_monomial_cliques(2*t,A)[1:end-1]
    # L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
    oddone = [[Y[[k, a[mc[k]]]] for a in moms₂ₜ[k]] for k in 1:p]
    # L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]   
    ondiag = [[ sqrt(A[i,i])*Y[[k, a[mc[k]]]] - Y[[k, (2*mom.eᵢ(n,i) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for i ∈ mc[k]]
    # L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
    offdiag = [[ A[i,j]*Y[[k, a[mc[k]]]] - Y[[k, (mom.eᵢ(n,i) + mom.eᵢ(n,j) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for i ∈ mc[k] for j ∈ setdiff(Set(mc[k]),Set([i]))]
    return [oddone..., offdiag..., ondiag...]
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
    return [[ Y[[k,(a + b + mom.eᵢ(n,i) + mom.eᵢ(n,j))[mc[k]]]] for a in momsₜ₋₁[k], b ∈ momsₜ₋₁[k]] for k ∈ 1:p for i ∈ mc[k] for j ∈ setdiff(Set(mc[k]),Set([i]))  ]
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

############################################### Utility #########################################
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Y, index_array) = map(α -> Y[α], index_array)
α_to_Lxᵅ(Y, index_array, mck) =  map(γ -> Y[[k,γ[mck]]], index_array)

function run_tests()
    @testset "sparse" begin
        # @test var_inds = get_var_inds(mc,t) 
        # @test sort(unique([ b  for (a,b) ∈ var_inds])) ==  sort(mom.make_mon_expo(length(mc[k]),2*t))
        # @test length(spar_pos_cons) == length(cat(get_maximal_cliques(M)...,dims=1))
        # @test length(spar_loc_cons) == length(cat([nz(m) for m in get_maximal_cliques(M)]...,dims=1))
    end
end

end



