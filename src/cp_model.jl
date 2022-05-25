module cp_model
using LinearAlgebra ; const la = LinearAlgebra
using JuMP
using MosekTools # The solver that we use.

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")


using .moments ; const mom = moments
using .constraints ; const con = constraints

export modelξₜᶜᵖⁱᵈ,
       modelξₜᶜᵖˢᵖ,
       computeξₜᶜᵖ,
       extract_moments


"""
    min Lᵧ(1)
    s.t Lᵧ(xxᵀ) = A
        Lᵧ([x][x]ᵀ) ⪰ 0
        Lᵧ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]
        Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0
        Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0       # Lᵧ(xᵢxⱼ[x]₂ₜ₋₂) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
        Lᵧ((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0
"""
function modelξₜᶜᵖⁱᵈ(A,t;G_con=true)
    n = size(A)[1]
    model = Model()
    @variable(model, Y[mom.make_mon_expo(n,2*t)] ) # Define variables in the moment matrix.
    # Lᵧ([x][x]ᵀ) ⪰ 0
        mom_matₜ_expo = mom.make_mon_expo((t,t),A)
        mom_matₜ      = con.α_to_Lxᵅ(Y, mom_matₜ_expo)
        @constraint(model, la.Symmetric(mom_matₜ) in PSDCone())
    # Lᵧ(xxᵀ) = A
        mom_mat₌₁ = mom.make_mon_expo(n,(1,1),isle=false)
        Lx_mom_mat₌₁ = con.α_to_Lxᵅ(Y, mom_mat₌₁)
        for  i in 1:n, j in i:n
            fix(Lx_mom_mat₌₁[i,j], A[i,j])
        end  
    # L(u) ≧ 0 for u ∈ [x]₂ₜ
    # L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
    # # L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
    #     for c in con.make_dag_con(A,t,Y)
    #         @constraint(model, c .≥ 0)
    #     end
    # # L(xᵢxⱼu) ⪰ 0 for for u ∈ [x]₂ₜ₋₂ and i,j ∈ [n] 
    #     for c in con.make_xx_con(A,t,Y)
    #         @constraint(model, c in PSDCone())
    #     end
    # Lᵧ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]
    # Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0
        for c in con.make_loc_con(A,t,Y)
            @constraint(model, c in PSDCone())
        end  
    # Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
        for c in con.make_ideal_con(A,t,Y)
            @constraint(model, c .== 0)
        end     
    # Lᵧ((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0
        if G_con
            println("----------------G-constraints are active")
            G_con = con.make_G_con(A,t,Y)
            @constraint(model, la.Symmetric(G_con) in PSDCone())
        end
    # Lᵧ(1)
        @objective(model, Min, mom_matₜ[1,1])
    return model
end

"""
    min ∑ₖ Lᵧₖ(1)
    s.t ∑ₖ Lᵧₖ(xxᵀ) = A , 
        Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]
        Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i ∈ Vₖ, k ∈ [p]
        Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]
        Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]
"""
function modelξₜᶜᵖˢᵖ(M,t;G_con=true,isWeak=false)
    model = Model()
    mc = get_maximal_cliques(M)
    spar_inds = con.make_spar_inds(mc,t)
    @variable(model, Y[spar_inds] ) 
    # Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]
        for m ∈ con.make_spar_mom_mat_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end
    # ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]
        for (c,v) ∈ con.make_spar_sec_ord_mom_con(M,Y)
            @constraint(model, c == v)
        end  
    # L(u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ for  k ∈ 1:p 
    # L((√Aᵢᵢ xᵢ - xᵢ²)⋅u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂, k ∈ 1:p, i ∈ mc[k]
    # # L((Aᵢⱼ  - xᵢxⱼ)⋅u) ≥ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂ for k ∈ 1:p  if i,j ∈ mc[k] 
    # for c in con.make_spar_dag_con(M,t,Y)
    #     @constraint(model, c .≥ 0)
    # end
    # # L(xᵢxⱼu) ⪰ 0 for u ∈ [x]ᵥₖ_₂ₜ₋₂ ; i,j ∈ mc[k], k ∈ [p] 
    # for c in con.make_spar_xx_con(M,t,Y)
    #     @constraint(model, c in PSDCone())
    # end
    # Lᵧₖ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i ∈ Vₖ, k ∈ [p]    
    # Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]
        for m ∈ con.make_spar_loc_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end    
    # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]
        if G_con
            println("----------------G-constraints are active")
            for m ∈ con.make_spar_G_con(M,t,Y,isWeak)
                @constraint(model, m in PSDCone())
            end    
        end  
    # ∑ Lᵧₖ(1)
        @objective(model, Min, con.make_obj(mc,Y))
    return model
end

"""This is where the ξₜᶜᵖ is calculated for matrix A """
function computeξₜᶜᵖ(model)
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model) 

    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end

""" Returns L(xᵅ) for α ∈ ℕⁿ₂ₜ """
function extract_moments(model)
    Y = model[:Y]
    indices = model[:Y].axes[1]    # {a ∈ model ⊆  ℕⁿ₂ₜ or  ℕᵛ₂ₜ}                                     
    values = JuMP.value.(Y).data    # {yₐ ∈ model}
    return hcat(indices, values)
end





end
