module nn_model
using LinearAlgebra ; const la = LinearAlgebra
using JuMP
using MosekTools # The solver that we use.

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")

using .moments ; const mom = moments


export modelξₜⁿⁿⁱᵈ,
       modelξₜⁿⁿˢᵖ,
       computeξₜⁿⁿ,
       get_ξₜⁿⁿ

"""
    min L(1)
    s.t L(xᵢxⱼ) = Mᵢⱼ..................................................(sec)
        L([x]ₜ[x]ₜᵀ) ⪰ 0................................................(mom)
        L((√Mₘₐₓxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [m + n].................(loc)
        L((Mᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j ∈ Eₘ .......................(loc)
        L(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j ∈ Eₘ .............................(xx)
        L([x]₂ₜ) ≥ 0 ...................................................(dag)
        L((√Mₘₐₓxᵢ-xᵢ²)[x]₂ₜ₋₂) ≥ 0 , i ∈ [m + n] ......................(dag)
        L((Mᵢⱼ-xᵢxⱼ)[x]₂ₜ₋₂) ≥ 0 ,i≠j ∈ Eₘ .............................(dag)
        L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all {i,j} s.t. Mᵢⱼ = 0...................(ideal)
"""
function modelξₜⁿⁿⁱᵈ(M,t,dag=true, ddag=false, xx = false) 
    m,n = size(M)
    model = Model()
    @variable(model, Y[mom.make_mon_expo(m+n,2*t)] ) # Define variables in the moment matrix.
    # L([x]ₜ[x]ₜᵀ) ⪰ 0
        Lxx = get_Lxx(M,t,Y)
        @constraint(model, Lxx in PSDCone())
    # L(xᵢxⱼ) = Mᵢⱼ
        for c in get_sec_ord_cons(M,Y)
            # fix(c[1], c[2])
            @constraint(model, c[1] == c[2])
        end
    # L((√Mₘₐₓxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [m + n]
        for c ∈ get_loc1_cons(M,t,Y)
            if t == 1
                @constraint(model, c[1] ≥ 0)
            else
                @constraint(model, c in PSDCone())
            end
        end 
    # L((Mᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j ∈ Eₘ 
        for c ∈ get_loc2_cons(M,t,Y)
            if t == 1
                @constraint(model, c[1] ≥ 0)
            else
                @constraint(model, c in PSDCone())
            end
        end     
    # L(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ E s.t. Mᵢⱼ = 0 
    for c ∈ get_ideal_cons(M,t,Y)
        # fix(c, 0.0)
        @constraint(model, c .== 0) 
    end
    # L(u) ≥ 0 for u ∈ [x]₂ₜ
    # L((√Mₖₖ xₖ - xₖ²)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂ 
    # L((Mₖₕ  - xₖxₕ)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂, k,h ∈ U x W  
    if dag  
        for c in get_dagger_cons(M,t,Y,ddag)
            @constraint(model, c .≥ 0)
        end
    end
    # L(1)
        @objective(model, Min, Lxx[1])
    return model
end
function get_Lxx(M,t,Y)
    mom_matₜ_expo = mom.make_mon_expo(mom.make_mon_expo(t,M,true))
    return a_to_yₐ(Y, mom_matₜ_expo)
end 
function get_sec_ord_cons(M,Y)
    m,n = size(M)
    mom_mat₌₁ = mom.make_mon_expo(m+n,(1,1),isle=false)
    Lx_mom_mat₌₁ = a_to_yₐ(Y, mom_mat₌₁)
    nze = mom.get_nonzero_entries(ones(m,n), false) # 
    return [(Lx_mom_mat₌₁[i,j+m], M[i,j] ) for (i,j) ∈ nze]
end
function get_loc1_cons(M,t,Y)
    m,n = size(M)
    sqrtMₘₐₓ = sqrt(maximum(M))
    # momₜ₋₁g(i) = mom.make_mon_expo(t-1,M,i,true) 
    momₜ₋₁g = mom.make_mon_expo(t-1,M,true) 
    Mₘₐₓ_terms(i) = [sqrtMₘₐₓ*Y[a+b+mom.eᵢ(m+n,i)] - Y[a+b+2*mom.eᵢ(m+n,i)] for a ∈ momₜ₋₁g, b ∈ momₜ₋₁g]
    return [Mₘₐₓ_terms(i) for i ∈ 1:(m+n)]
end
function get_loc2_cons(M,t,Y)
    m,n = size(M)
    momₜ₋₁ = mom.make_mon_expo(t-1, M, true)
    nze = mom.get_nonzero_entries(M, false)
    M_terms(i,j) = [M[i,j-m]*Y[a+b] - Y[a+b+mom.eᵢ(m+n,i,j)] for a ∈ momₜ₋₁, b ∈ momₜ₋₁]
    return [M_terms(i,j+m) for (i,j) ∈ nze]
end
function get_ideal_cons(M,t,Y)
    m,n = size(M)
    mom₂ₜ₋₂ = mom.make_mon_expo(2t-2, M, true) 
    ze = mom.get_zero_entries(M, false)
    return [Y[a+mom.eᵢ(m+n,i,j+m)] for (i,j) ∈ ze, a ∈ mom₂ₜ₋₂]
end
function get_dagger_cons(M,t,Y,dd=false)
    m,n = size(M)
    sqrtMₘₐₓ = sqrt(maximum(M))
    mom₂ₜ   = mom.make_mon_expo(2t, M, true)
    mom₂ₜ₋₂ = mom.make_mon_expo(2t-2, M, true)
    nze = mom.get_nonzero_entries(ones(m,n), false)
    # L(u) ≥ 0 for u ∈ [x]₂ₜ
    oddone = [[ Y[a] for a in mom₂ₜ]]
    # L((√Mₖₖ xₖ - xₖ²)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂ 
    ondiag = [[ sqrtMₘₐₓ*Y[mom.eᵢ(m+n,i) + a] - Y[2*mom.eᵢ(m+n,i) + a]  for a in mom₂ₜ₋₂] for i ∈ 1:m+n ]
    # L((Mₖₕ  - xₖxₕ)⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂, k,h ∈ U x W 
    offdiag = [[ M[i,j]*Y[a] - Y[mom.eᵢ(m+n,i,j+m) + a]  for a in mom₂ₜ₋₂] for (i,j) ∈ nze]
    return   dd ? [oddone..., ondiag...,offdiag...] : offdiag  # , 
end
function get_xx_cons(M,t,Y)
    m,n = size(M)
    mom₂ₜ₋₂ = mom.make_mon_expo(2t-2, M, true)
    nze = mom.get_nonzero_entries(ones(m,n), false)
    # L(xₖxₕ⋅u) ≥ 0 for u ∈ [x]₂ₜ₋₂, k,h ∈ U x W 
    return   [[Y[mom.eᵢ(m+n,i,j+m) + a]  for a in mom₂ₜ₋₂] for (i,j) ∈ nze]
end

########################################################### Sparse ################################################################

"""
    min ∑ₖ Lₖ(1)
    s.t ∑ₖ Lₖ(xᵢxⱼ) = Mᵢⱼ................................................(sec)
        Lₖ ([x(Vₖ)]ₜ[x(Vₖ)]ₜᵀ) ⪰ 0 ∀ k ∈ [p]..............................(mom)
        Lₖ ((√Mₘₐₓxᵢ-xᵢ²)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0 , i ∈ Vₖ, ∀ k ∈ [p]...(loc)
        Lₖ((Mᵢⱼ-xᵢxⱼ)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0 ,{i,j} ⊆ Vₖ, ∀ k ∈ [p] ...(loc)
        L([x]₂ₜ) ≥ 0 ....................................................(dag)
        L((√Mₘₐₓxᵢ-xᵢ²)[x]₂ₜ₋₂) ≥ 0 , i ∈ [m + n] .......................(dag)
        L((Mᵢⱼ-xᵢxⱼ)[x]₂ₜ₋₂) ≥ 0 ,i≠j ∈ Eₘ ..............................(dag)
        Lₖ(xᵢxⱼ[x(Vₖ)]₂ₜ₋₂) ≥ 0 for all i,j ∈ Vₖ s.t. Mᵢⱼ ≠ 0.............(xx)
"""
function modelξₜⁿⁿˢᵖ(M,t,dag=true, ddag=false, xx=false)
    mc = mom.get_maximal_cliques(M,true)
    model = Model()
    spar_inds = make_spar_inds(mc,t)
    @variable(model, Y[spar_inds] ) 
        for c ∈ get_L_kxxV_k(M,Y,t)
            @constraint(model, c in PSDCone())
        end
    # ∑ₖ Lₖ(xᵢxⱼ) = Mᵢⱼ 
        for c in get_sparse_sec_ord_cons(M,Y)
            @constraint(model, c[1] == c[2])
        end
    # Lₖ((√Mₘₐₓxᵢ-xᵢ²)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0 for i ∈ Vₖ for k ∈ [p]
        for c ∈ get_sparse_loc1_cons(M,t,Y)
            if t == 1
                @constraint(model, c[1] ≥ 0)
            else
                @constraint(model, c in PSDCone())
            end
        end
    # Lₖ((Mᵢⱼ-xᵢxⱼ)[x(Vₖ)]ₜ₋₁[x(Vₖ)]ₜ₋₁ᵀ) ⪰ 0 ,{i,j} ⊆ Vₖ, ∀ k ∈ [p]
        if t == 1
            for c in get_sparse_loc2_cons(M,t,Y)
                @constraint(model, c[1] ≥ 0)
            end
        else
            for c in get_sparse_loc2_cons(M,t,Y)
                @constraint(model, c in PSDCone())
            end
        end
    # Lₖ(u) ≥ 0 for u ∈ [x(Uₖ∪Wₖ)]₂ₜ for k ∈ [p]
    # Lₖ((√Mᵢᵢ xᵢ - xᵢ²)⋅u) ≥ 0 for i ∈ Uₖ∪Wₖ,  u ∈ [x(Uₖ∪Wₖ)]₂ₜ₋₂, k ∈ [p]
    # L((Mᵢⱼ  - xᵢxⱼ)⋅u) ≥ 0 for i,j ∈ Uₖ × Wₖ, u ∈ [x(Uₖ∪Wₖ)]₂ₜ₋₂, k ∈ [p] 
    if dag 
        for c in get_sparse_dagger_cons(M,t,Y,ddag)
            @constraint(model, c .≥ 0)
        end
    end
    #    
    @objective(model, Min, sum([Y[[k,zeros(Int,length(mc[k]))]] for k in 1:length(mc)]))
    return model
end
make_spar_inds(mc,t) = [ [k,m] for  k ∈ 1:length(mc) for m in mom.make_mon_expo(length(mc[k]), 2*t) ]
function get_L_kxxV_k(M,Y,t)
    mc = mom.get_maximal_cliques(M,true)
    p = length(mc)
    Iᵏs = mom.get_monomial_cliques((t,t),M,true)[1:end-1] # Lₖ ([x(Vₖ)]ₜ[x(Vₖ)]ₜᵀ) ⪰ 0 ∀ k ∈ [p]
    Iᵏs_cut = [map(a -> [k, a[mc[k]] ], Iᵏs[k]) for k ∈ 1:p]
    return [a_to_yₐ(Y, Iᵏ) for Iᵏ ∈ Iᵏs_cut]
end
function get_sparse_sec_ord_cons(M,Y)
    m,n  = size(M)
    mc   = mom.get_maximal_cliques(M,true)
    p    = length(mc)
    xx₁  = mom.make_mon_expo(m+n,(1,1),isle=false)[1:m, (m+1):(m+n)]   # [x]₌₁[x]₌₁ᵀ
    Iᵏs₁ = mom.get_monomial_cliques((1,1), M,true)[1:end-1]
    ∑Lₖxx₁ = sum([map(a -> (a ∈ Iᵏs₁[k]) ? Y[[k,a[mc[k]]]] : 0, xx₁) for k ∈ 1:p])
    return [(∑Lₖxx₁[i,j], M[i,j])  for  i in 1:m, j in 1:n if M[i,j] != 0.0 ]
end
function get_sparse_loc1_cons(M,t,Y)
    m,n = size(M)
    sqrtMₘₐₓ = sqrt(maximum(M))
    mc = mom.get_maximal_cliques(M, true)
    p = length(mc)
    Iᵏsₜ₋₁ = mom.get_monomial_cliques(t-1, M, true)[1:end-1]  # [x(Vₖ)]ₜ₋₁ for k ∈ [p]        
    xVₖ = [map(a->a[mc[k]],Iᵏsₜ₋₁[k])  for k ∈ 1:p ]          # [x(Vₖ)]ₜ₋₁ with exponents trimmed to be in Vₖ for k ∈ [p]
    eᵢVₖ(k,i)     = mom.eᵢ(m+n,i)[mc[k]]
    lht(a,b,k,i)  = sqrtMₘₐₓ*Y[[k, a+b+eᵢVₖ(k,i)]]         # Lₖ (√Mₘₐₓxᵢxᵅ⁺ᵝ)
    rht(a,b,k,i)  = Y[[k, a+b+2*eᵢVₖ(k,i)]]                # Lₖ (xᵢ²xᵅ⁺ᵝ) 
    return [[lht(a,b,k,i) - rht(a,b,k,i) for a ∈ xVₖ[k], b ∈ xVₖ[k]] for k ∈ 1:p  for i ∈ mc[k]]
end
function get_sparse_loc2_cons(M,t,Y)
    m,n = size(M)
    mc = mom.get_maximal_cliques(M, true)
    p = length(mc)
    Iᵏsₜ₋₁ = mom.get_monomial_cliques(t-1, M, true)[1:end-1]  # [x(Vₖ)]ₜ₋₁ for k ∈ [p]        
    xVₖ = [map(a->a[mc[k]],Iᵏsₜ₋₁[k])  for k ∈ 1:p ] # [x(Vₖ)]ₜ₋₁ with exponents trimmed to be in Vₖ for k ∈ [p]

    eᵢVₖ(k,i)    = mom.eᵢ(m+n,i)[mc[k]]
    eᵢⱼVₖ(k,i,j) = mom.eᵢ(m+n,i,j)[mc[k]]

    lht(a,b,k,i,j)  = M[i,j-m]*Y[[k, a+b]]              # Lₖ (Mᵢⱼxᵅ⁺ᵝ)
    rht(a,b,k,i,j)  =          Y[[k, a+b+eᵢⱼVₖ(k,i,j)]] # Lₖ (xᵢxⱼxᵅ⁺ᵝ)

    return [[lht(a,b,k,i,j) - rht(a,b,k,i,j) for a ∈ xVₖ[k], b ∈ xVₖ[k]]  for k ∈ 1:p for (i,j) ∈ mom.get_edges(mc[k],M) ]
end
function get_sparse_dagger_cons(M,t,Y,dd=false)
    m,n = size(M)
    sqrtMₘₐₓ = sqrt(maximum(M))
    mom₂ₜ   = mom.get_monomial_cliques(2t,M,true)[1:end-1]
    mom₂ₜ₋₂ = mom.get_monomial_cliques(2t-2,M,true)[1:end-1]
    mc = mom.get_maximal_cliques(M,true)
    p = length(mc)
    # Lₖ(u) ≥ 0 for u ∈ [x(Uₖ∪Wₖ)]₂ₜ for k ∈ [p]
    oddone =  [ [Y[[k,a[mc[k]]]]  for a in mom₂ₜ[k]] for k ∈ 1:p ] 
    # Lₖ((√Mᵢᵢ xᵢ - xᵢ²)⋅u) ≥ 0 for i ∈ Uₖ∪Wₖ,  u ∈ [x(Uₖ∪Wₖ)]₂ₜ₋₂, k ∈ [p]
    ondiag =  [ [sqrtMₘₐₓ*Y[[k,(mom.eᵢ(m+n,i) + a)[mc[k]]]] - Y[[k,(2*mom.eᵢ(m+n,i) + a)[mc[k]]]]  for a in mom₂ₜ₋₂[k]] for k ∈ 1:p for i ∈ mc[k]]
    # L((Mᵢⱼ  - xᵢxⱼ)⋅u) ≥ 0 for i,j ∈ Uₖ × Wₖ, u ∈ [x(Uₖ∪Wₖ)]₂ₜ₋₂, k ∈ [p] 
    offdiag = [ [M[i,j-m]*Y[[k,a[mc[k]]]] - Y[[k,(mom.eᵢ(m+n,i,j)+a)[mc[k]]]]  for a in mom₂ₜ₋₂[k] ] for k ∈ 1:p for (i,j) ∈ mom.get_edges(mc[k],M)]
    return  dd ? [oddone..., ondiag..., offdiag...] : offdiag  # ,  
end
function get_sparse_xx_cons(M,t,Y)
    error() 
 end

"""This is where the ξₜᶜᵖ is calculated for matrix A """
function computeξₜⁿⁿ(model)
    set_optimizer(model, Mosek.Optimizer)
    # set_silent(model)
    optimize!(model) 

    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end

function get_ξₜⁿⁿ(M,t,flavour)
    @assert !contains(flavour,"iddag")
    contains(flavour,"dag")  ? dag   =true : dag  =false
    contains(flavour,"ddag") ? ddag  =true : ddag =false 
    contains(flavour,"xx") ? xx = true : xx = false
    if      contains(flavour,"id")
        mod = modelξₜⁿⁿⁱᵈ(M, t, dag, ddag, xx)
    elseif  contains(flavour,"sp")
        mod = modelξₜⁿⁿˢᵖ(M, t, dag, ddag, xx)
    else
        error("Incorrect model specification")
    end
    ξₜⁿⁿ = computeξₜⁿⁿ(mod) ;
    ex_mom = nothing #extract_moments(ξₜᶜᵖ)
    return ξₜⁿⁿ, ex_mom 
end

############################################### Utility Utility #########################################
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
a_to_yₐ(Y, i_arr)    = map(a -> Y[a], i_arr)
a_to_yₐ(Y, i_arr, c) = map(a -> Y[[k,a[c]]], i_arr)

end





