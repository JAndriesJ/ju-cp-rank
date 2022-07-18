module cp_model
using LinearAlgebra ; const la = LinearAlgebra
using JuMP
using MosekTools # The solver that we use.

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
using .moments ; const mom = moments

export get_ξₜᶜᵖ,  
       extract_moments

################################################################## Dense #######################################################################

"""
    min L(1)
    s.t L(xxᵀ) = A
        L([x][x]ᵀ) ⪰ 0
        L((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]
        L((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0
        L(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0       
        L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
        L((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0
"""
function modelξₜᶜᵖⁱᵈ(A,t;G_con=true,dag=true,ddag=false,xx=false)
    n = size(A)[1]
    model = Model()
    # L([x]₂ₜ)------------------------------------------------------------------------------------------------------------------------
    @variable(model, Y[mom.make_mon_expo(n, 2*t)] ) 
    #L([x][x]ᵀ) ⪰ 0------------------------------------------------------------------------------------------------------------------
    @constraint(model, make_mom_mat_con(A,t,Y) in PSDCone())  
    #L(xxᵀ) = A------------------------------------------------------------------------------------------------------------------------
    for (c,v) ∈ make_sec_ord_mom_con(A,Y) @constraint(model, c == v) end 
    #L(u) ≧ 0 for u ∈ [x]₂ₜ 
    # L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂ 
    # L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂ ---------------------------------------------------------------------------------------
    if dag if t != 1  
        println("----------------dag-constraints are active") 
        for c in make_dag_con(A,t,Y,ddag) @constraint(model, c .≥ 0) end
    end end
    # L(xᵢxⱼu) ⪰ 0 for for u ∈ [x]₂ₜ₋₂ and i,j ∈ [n]-----------------------------------------------------------------------------------
    if xx if t != 1 
            println("----------------xx-constraints are active")
            for c in make_xx_con(A,t,Y)  @constraint(model, c in PSDCone())  end  
    end end
    #L((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]----L((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0-------------------------------------------
    for c in make_loc_con(A,t,Y)
        t == 1 ? @constraint(model, c[1] ≥ 0 ) : @constraint(model, c in PSDCone())
    end  
    #L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0----------------------------------------------------------------------------------
    for c in make_ideal_con(A,t,Y) # 
        t == 1 ? fix(c[1], 0.0 ) : @constraint(model, c .== 0) 
    end
    #L((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0-----------------------------------------------------------------------------------------------------
    if G_con 
        println("----------------G-constraints are active")
        @constraint(model, make_G_con(A,t,Y) in PSDCone())
    end
    #min L(1)--------------------------------------------------------------------------------------------------------------------------   
    @objective(model, Min, Y[zeros(Int,n)])  
    return model
end

make_mom_mat_con(A,t,Y) =  a_to_yₐ(Y,mom.make_mon_expo(mom.make_mon_expo(t,A)))

function make_sec_ord_mom_con(A,Y)
    Lx_mom_mat₌₁ = a_to_yₐ(Y, mom.make_mon_expo(size(A)[1], (1,1), isle=false))
    nze = mom.get_nonzero_entries(ones(size(A)...))
    return [(Lx_mom_mat₌₁[ij...], A[ij...]) for ij ∈ nze]
end

"""(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ)  ;  (M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function make_loc_con(A,t,Y)
    n = size(A)[1]
    nze = mom.get_nonzero_entries(A)
    momₜ₋₁ = mom.make_mon_expo(t-1, A)
    momₜ₋₁g(j) = mom.make_mon_expo(t-1,A,j)
    if_diag(k) = [sqrt(A[k,k])*Y[α+β+mom.eᵢ(n,k)] - Y[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁g(k), β ∈ momₜ₋₁g(k)]
    if_off_diag(k,h) = [A[k,h]*Y[α+β] - Y[α+β+mom.eᵢ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    return [e[1]==e[2] ? if_diag(e[1]) : if_off_diag(e[1],e[2]) for e in nze]
end

"""M(G ⊗ L) ⪰ 0 constraints input: A(data matrix),t(Integer),Y(JuMP variable) Assumption: G = A-[x]₌₁[x]₌₁ᵀ output: A⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ) """
function make_G_con(A,t,Y)
    n = size(A)[1]
    LMBexp₁    = mom.make_mon_expo(n,(1,1),isle=false)          # exponents of [x]₌₁[x]₌₁ᵀ
    LMBexpₜ₋₁   = mom.make_mon_expo(mom.make_mon_expo(t-1,A))    # exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁      = a_to_yₐ(Y,LMBexpₜ₋₁)                           # L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp₁ₜ₋₁  = mom.expo_kron(LMBexp₁,LMBexpₜ₋₁)                # exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB₁ₜ₋₁     = a_to_yₐ(Y,LMBexp₁ₜ₋₁)                          # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB₁ₜ₋₁                              # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    return G_con
end

"""L(xᵢxⱼ[x]₂ₜ₋₂) = 0  for all {i,j} s.t. Mᵢⱼ = 0"""
function make_ideal_con(A,t,Y) 
    n = size(A)[1] 
    mom₂ₜ₋₂ = mom.make_mon_expo(2t-2, A)
    nze = mom.get_zero_entries(A)
    eᵢⱼs = map( e -> mom.eᵢ(n,e[1],e[2]), nze)
    return map(eᵢⱼ -> [Y[a+eᵢⱼ] for a ∈ mom₂ₜ₋₂ ], eᵢⱼs)
end

"""input: A(Data matrix),t(integer),Y(JuMP variable)   output: L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)"""
function make_dag_con(A,t,Y, dd=false)
    n = size(A)[1]
    mom₂ₜ   = mom.make_mon_expo(2*t, A)  # [x]₂ₜ
    mom₂ₜ₋₂ = mom.make_mon_expo(2*t-2, A) # [x]₂ₜ₋₂
    oddone = [[ 1*Y[α] for α in mom₂ₜ]] # L(u) ≥ 0 for u ∈ [x]₂ₜ
    ondiag = [[ sqrt(A[i,i])*Y[mom.eᵢ(n,i) + a] - Y[2*mom.eᵢ(n,i) + a]  for a in mom₂ₜ₋₂] for i ∈ 1:n ] # L((√Aₖₖ xₖ - xₖ²)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂
    offdiag = [[ A[i,j]*Y[a] - Y[mom.eᵢ(n,i) + mom.eᵢ(n,j) + a]  for a in mom₂ₜ₋₂] for i ∈ 1:n for j ∈ i+1:n if A[i,j] != 0] # L((Aᵢⱼ  - xᵢxⱼ)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂, {i,j} ∈ E
    dd ?  println("----------------double dag-constraints are active") : 0
    return  dd ? [oddone...,offdiag...,  ondiag...] : offdiag 
end
   
""" input: A(Data matrix), t(integer), Y(JuMP variable)   output: L(xᵢxⱼu) ⪰ 0 for u ∈ [x]₂ₜ₋₂ and i,j ∈ [n] """
function make_xx_con(A,t,Y)
    n = size(A)[1]
    mom₁ₜ₋₁ = mom.make_mon_expo(t-1,A)
    return [[ Y[a + b + mom.eᵢ(n,i) + mom.eᵢ(n,j)]  for a in mom₁ₜ₋₁, b ∈ mom₁ₜ₋₁] for i ∈ 1:n for j ∈ i+1:n if A[i,j] != 0.0] # L(xᵢxⱼu) for g ∈ {1} ∪ Sᶜᵖ_A and i,j ∈ [n]
end

################################################################## Sparse #######################################################################

"""
    min ∑ₖ Lᵧₖ(1)
    s.t ∑ₖ Lᵧₖ(xxᵀ) = A , 
        Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]
        Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i ∈ Vₖ, k ∈ [p]
        Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]
        Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]
"""
function modelξₜᶜᵖˢᵖ(M,t;G_con=true,dag=false,ddag=false,xx=false,isWeak=false)
    model = Model()
    mc = get_maximal_cliques(M)
    spar_inds = make_spar_inds(mc,t)
    #L([x(Vₖ)]₂ₜ)-------------------------------------------------------------------------------------------------------------------------   
    @variable(model, Y[spar_inds] ) 
    #Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]----------------------------------------------------------------------------------------------------------
    for m ∈ make_spar_mom_mat_con(M,t,Y) @constraint(model, m in PSDCone()) end
    #∑ Lᵧₖ(xxᵀ) = A , k ∈ [p] -----------------------------------------------------------------------------------------------------------      
    for (c,v) ∈ make_spar_sec_ord_mom_con(M,Y) @constraint(model, c == v) end 
    #L(u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p----L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ Vₖ--------------------------------
    #L((Aᵢⱼ  - xᵢxⱼ)⋅u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ Vₖ ------------------------------------------------------------------
    if dag if t != 1
            println("----------------dag-constraints are active")
            for c in make_spar_dag_con(M,t,Y,ddag) @constraint(model, c .≥ 0) end
    end end
    #Lₖ(xᵢxⱼu) ⪰ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂ ; i,j ∈ mc[k], k ∈ [p] ------------------------------------------------------------------------------
    if xx if t != 1 
            println("----------------xx-constraints are active")
            for c in make_spar_xx_con(M,t,Y) @constraint(model, c in PSDCone()) end
    end end
    #Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i ∈ Vₖ, k ∈ [p]----Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]----------------------
    for c ∈ make_spar_loc_con(M,t,Y)
        t == 1 ? @constraint(model, c[1] ≥ 0 ) : @constraint(model, c in PSDCone())
    end  
    #Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]-------------------------------------------------------------------------------------         
    if G_con
        println("----------------G-constraints are active")
        for m ∈ make_spar_G_con(M,t,Y,isWeak)  @constraint(model, m in PSDCone()) end    
    end
    #∑ Lᵧₖ(1) ------------------------------------------------------------------------------------------------------------------------------         
        @objective(model, Min, make_obj(mc,Y))
    return model
end

make_spar_inds(mc,t) = [ [k,m] for  k ∈ 1:length(mc) for m in mom.make_mon_expo(length(mc[k]), 2*t) ]

""" ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]"""
function make_spar_sec_ord_mom_con(M,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    xx₁ = mom.make_mon_expo(n,(1,1),isle=false)
    Iᵏs₁ = mom.get_monomial_cliques((1,1), M)[1:end-1]
    K = sum([map(α -> (α ∈ Iᵏs₁[k]) ? Y[[k,α[mc[k]]]] : 0, xx₁) for k ∈ 1:p])
    nze = mom.get_nonzero_entries(M)
    return [(K[ij...], M[ij...])  for ij ∈ nze]
end

"""Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]""" 
function make_spar_mom_mat_con(M,t,Y)
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = mom.get_monomial_cliques((t,t),M)[1:end-1]
    Iᵏs_cut = [map(a -> [k, a[mc[k]] ], Iᵏs[k]) for k ∈ 1:p]
    return [a_to_yₐ(Y, Iᵏ) for Iᵏ ∈ Iᵏs_cut]
end

""" L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
    L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]
    L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
"""
function make_spar_dag_con(A,t,Y, dd=false)
    n = size(A)[1]
    mc = moments.get_maximal_cliques(A)
    p = length(mc)
    moms₂ₜ₋₂ = moments.get_monomial_cliques(2*t-2,A)[1:end-1]
    moms₂ₜ   = moments.get_monomial_cliques(2*t,A)[1:end-1]

    oddone = [[Y[[k, a[mc[k]]]] for a in moms₂ₜ[k]] for k in 1:p] # L(u) for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
    ondiag = [[sqrt(A[i,i])*Y[[k, (mom.eᵢ(n,i)+a)[mc[k]]]] - Y[[k, (2*mom.eᵢ(n,i) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for i ∈ mc[k]] # L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k] 
    offdiag = [[A[i,j]*Y[[k, a[mc[k]]]] - Y[[k, (mom.eᵢ(n,i) + mom.eᵢ(n,j) + a)[mc[k]]]]  for a in moms₂ₜ₋₂[k]] for k in 1:p for (i,j) ∈ mom.get_edges(mc[k],false)] # L((Aᵢⱼ  - xᵢxⱼ)⋅u) for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
    # warning("This could be wrong")                                                                                                             
    dd ?  println("----------------double dag-constraints are active") : 0
    return  dd ? [oddone...,offdiag...,  ondiag...] : offdiag 
end

"""L(xᵢxⱼu) for u ∈ [x]ᵥₖ_₂ₜ₋₂ ; i,j ∈ mc[k], k ∈ [p] """
function make_spar_xx_con(A,t,Y)
    n = size(A)[1]
    mc = moments.get_maximal_cliques(A)
    p = length(mc)
    momsₜ₋₁ = moments.get_monomial_cliques(t-1,A)
    return [
            [ Y[[k,(a + b + mom.eᵢ(n,i) + mom.eᵢ(n,j))[mc[k]]]]
               for a in momsₜ₋₁[k], b ∈ momsₜ₋₁[k] 
                                                 ] for k ∈ 1:p for (i,j) ∈ mom.get_edges(mc[k],false)  ]
end

""" Lₖ((√Aᵢᵢxᵢ-xᵢ²)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ)  ⪰ 0, i ∈ Vₖ  ;  Lₖ((Aᵢⱼ-xᵢxⱼ)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0, i,j ∈ Vₖ, k ∈ [p] """
function make_spar_loc_con(M,t,Y)
    n = size(M)[1]
    mc = get_maximal_cliques(M)
    p = length(mc)
    Iᵏs = get_monomial_cliques(t-1,M)

    term(e,k,a,b)     = e[1]==e[2] ? d_term(k,e,a,b) : od_term(k,e,a,b)
    d_term(k,e,a,b)   = sqrt(M[e...])*Y[ [k, (a+b+eᵢ(n,e[1]) )[mc[k]]] ] - Y[ [k, (a+b + eᵢ(n,e[1],e[2]))[mc[k]]]] # Lₖ((√Mᵢᵢxᵢ-xᵢ²)xᵅ)
    od_term(k,e,a,β)  =      M[e...] *Y[ [k, (a+β            )[mc[k]]] ] - Y[ [k, (a+β + eᵢ(n,e[1],e[2]))[mc[k]]]] # Lₖ((Aᵢⱼ-xᵢxⱼ)xᵅ)

    arr = [[[term(e,k,a,b)  for a ∈ Iᵏs[k], b ∈ Iᵏs[k]]              # iterate over monomial exponent ass. to a clique
                            for e ∈ moments.get_edges(sort(mc[k]))]  # iterate over each pair of element of the clique
                            for k ∈ 1:p]                             # iterate over each clique
    return cat(arr...,dims=1) 
end

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

################################################################## Utility #######################################################################
"""This is where the ξₜᶜᵖ is calculated for matrix A """
function computeξₜᶜᵖ(model)
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model) 

    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end

"""
Computes one of the many flavours of ξₜᶜᵖ
"""
function get_ξₜᶜᵖ(M, t, flavour) 
    @assert !contains(flavour,"iddag") 
    contains(flavour,"dag")  ? dag   =true : dag  =false
    contains(flavour,"ddag") ? ddag  =true : ddag  =false 
    contains(flavour,"xx")   ? xx    =true : xx   =false
    contains(flavour,"G")    ? G_con =true : G_con=false
    if      contains(flavour,"id")
        mod = modelξₜᶜᵖⁱᵈ(M, t, G_con=G_con, dag=dag, ddag=ddag, xx=xx)
    elseif  contains(flavour,"wsp")
        mod = modelξₜᶜᵖˢᵖ(M, t, G_con=G_con, dag=dag, ddag=ddag, xx=xx, isWeak=true)
    elseif  contains(flavour,"sp")
        mod = modelξₜᶜᵖˢᵖ(M, t, G_con=G_con, dag=dag, ddag=ddag, xx=xx)
    else
        error("Incorrect model specification")
    end
    ξₜᶜᵖ = computeξₜᶜᵖ(mod)
    ex_mom = extract_moments(ξₜᶜᵖ)
    return ξₜᶜᵖ, ex_mom
end  


# function get_ξₜᶜᵖ(args, save_path::String)
#     ξₜᶜᵖ, mom, s = capture_solver_output(get_ξₜᶜᵖ, (args...)) # args = (M ,t, flavour)
#     write_solver_output(s, save_path)
#     return ξₜᶜᵖ, mom
# end

# function capture_solver_output(func,args)
#     original_stdout = stdout;
#     (rd, _) = redirect_stdout();
#     ξₜᶜᵖ, mom = func(args...)
#     s = []
#     for rl in eachline(rd)
#         push!(s,rl)
#         if contains(rl,"Objective:")
#             break
#         end
#     end
#     redirect_stdout(original_stdout);
#     return ξₜᶜᵖ, mom, s
# end
# function write_solver_output(s,save_path)
#     touch(save_path)
#     open(save_path, "w") do io
#         for line in s
#             write(io, line*"\n")
#         end
#     end;
# end

""" Returns L(xᵅ) for α ∈ ℕⁿ₂ₜ """
function extract_moments(model)
    Y = model[:Y]
    indices = model[:Y].axes[1]    # {a ∈ model ⊆  ℕⁿ₂ₜ or  ℕᵛ₂ₜ}                                     
    values = JuMP.value.(Y).data   # {yₐ ∈ model}
    return hcat(indices, values)
end

############################################### Utility Utility #########################################
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
a_to_yₐ(Y, i_arr)    = map(a -> Y[a], i_arr)
a_to_yₐ(Y, i_arr, c) = map(a -> Y[[k,a[c]]], i_arr)
end


# la.Symmetric(G_con)
# mom_matₜ_expo = mom.make_mon_expo(mom.make_mon_expo(t,A))
