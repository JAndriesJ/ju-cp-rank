module cp_model
using LinearAlgebra ; const la = LinearAlgebra
using JuMP
using MosekTools # The solver that we use.

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
include(proj_dir*"GMP_constraints.jl")


using .moments ; const mom = moments
using .constraints ; const con = constraints
using .GMP_constraints ; const gmpcon = GMP_constraints

export modelξₜᶜᵖ,
       GMP_dens_modelξₜᶜᵖ,
       GMP_spar_modelξₜᶜᵖ,
       computeξₜᶜᵖ

"""
    min Lᵧ(1)
    s.t Lᵧ(xxᵀ) = A
        Lᵧ([x][x]ᵀ) ⪰ 0
        Lᵧ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]
        Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0
        Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
        Lᵧ((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0
"""
function modelξₜᶜᵖ(A,t;G_con=true)
    n = size(A)[1]
    model = Model()
    @variable(model, Lx[mom.make_mon_expo(n,2*t)] ) # Define variables in the moment matrix.
    # Lᵧ([x][x]ᵀ) ⪰ 0
        mom_matₜ_expo = mom.make_mon_expo((t,t),A)
        mom_matₜ      = con.α_to_Lxᵅ(Lx, mom_matₜ_expo)
        @constraint(model, la.Symmetric(mom_matₜ) in PSDCone())
    # Lᵧ(xxᵀ) = A
        mom_mat₌₁ = mom.make_mon_expo(n,(1,1),isle=false)
        Lx_mom_mat₌₁ = con.α_to_Lxᵅ(Lx, mom_mat₌₁)
        for  i in 1:n, j in i:n
            fix(Lx_mom_mat₌₁[i,j], A[i,j])
        end     
    # Lᵧ((√Aᵢᵢxᵢ-xᵢ²)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , i ∈ [n]
    # Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 ,i≠j,Aᵢⱼ ≠ 0
        for c in con.make_loc_con(A,t,Lx)
            @constraint(model, c in PSDCone())
        end  
    # Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
        for c in con.make_ideal_con(A,t,Lx)
            @constraint(model, c .== 0)
        end     
    # Lᵧ((A-xxᵀ) ⊗ [x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0
        if G_con
            println("----------------G-constraints are active")
            G_con = con.make_G_con(A,t,Lx)
            @constraint(model, la.Symmetric(G_con) in PSDCone())
        end
    # Lᵧ(1)
        @objective(model, Min, mom_matₜ[1,1])
    return model
end

#################### GMP models ##########################
"""
    min Lᵧ(1)
    s.t Lᵧ(xxᵀ) = A
        Lᵧ([x][x]ᵀ) ⪰ 0
        Lᵧ(xᵢ [x]ₜ₋₁[x]ᵀₜ₋₁ ) for i∈ [n]
        Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , Aᵢⱼ ̸= 0, k ∈ [p]
        Lᵧ([x]ₜ₋₁[x]ₜ₋₁ᵀ ⊗ (A-xxᵀ)) ⪰ 0
        Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0
"""
function GMP_dens_modelξₜᶜᵖ(M,t;G_con=true)
    model = Model()
    dens_inds = gmpcon.make_dens_inds(M,t)
    @variable(model, Y[dens_inds] ) # Define variables in the moment matrix.
    dens_mom_mat = gmpcon.make_dens_mom_mat_con(M,t,Y)
    # Lᵧ([x][x]ᵀ) ⪰ 0
        @constraint(model, dens_mom_mat in PSDCone()) 
    # Lᵧ(xxᵀ) = A
        @constraint(model, gmpcon.make_dens_sec_ord_mom_con(M,Y) .== M)   
    # Lᵧ(xᵢ [x]ₜ₋₁[x]ᵀₜ₋₁ ) for i∈ [n]
        for m ∈ gmpcon.make_dens_pos_con(M,t,Y)
            @constraint(model, m in PSDCone()) 
        end     
    # Lᵧ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁[x]ₜ₋₁ᵀ) ⪰ 0 , Aᵢⱼ ̸= 0, k ∈ [p]
        for m ∈ gmpcon.make_dens_loc_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end  
    # Lᵧ(xᵢxⱼ[x]ₜ₋₁[x]ₜ₋₁ᵀ) = 0 for all i≠j ∈ [n] s.t. Mᵢⱼ = 0    
        for m ∈ gmpcon.make_dens_ideal_con(M,t,Y)
            @constraint(model, m  .== 0)
        end     
    # Lᵧ([x]ₜ₋₁[x]ₜ₋₁ᵀ ⊗ (A-xxᵀ)) ⪰ 0    
        if G_con
            println("----------------G-constraints are active")
            @constraint(model, gmpcon.make_dens_G_con(M,t,Y) in PSDCone())   
        end  
    # Lᵧ(1)
        @objective(model, Min, dens_mom_mat[1,1])
    return model
end


"""
    min ∑ Lᵧₖ(1)
    s.t ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]
        Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]
        Lᵧₖ(xᵢ [x]ₜ₋₁ᵥₖ [x]ᵀₜ₋₁ᵥₖ ) for i∈ Vₖ and k ∈ [p]
        Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]
        Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]
"""
function GMP_spar_modelξₜᶜᵖ(M,t;G_con=true)
    model = Model()
    mc = get_maximal_cliques(M)
    spar_inds = gmpcon.make_spar_inds(mc,t)
    @variable(model, Y[spar_inds] ) 
    # Lᵧₖ([x][x]ᵀ) ⪰ 0 , k ∈ [p]
        for m ∈ gmpcon.make_spar_mom_mat_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end
    # ∑ Lᵧₖ(xxᵀ) = A , k ∈ [p]
        for (c,v) ∈ gmpcon.make_spar_sec_ord_mom_con(M,Y)
            @constraint(model, c == v)
        end
       
    # Lᵧₖ(xᵢ [x]ₜ₋₁ᵥₖ [x]ᵀₜ₋₁ᵥₖ ) for i∈ Vₖ and k ∈ [p]
        for m ∈ gmpcon.make_spar_pos_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end     
    # Lᵧₖ((Aᵢⱼ-xᵢxⱼ)[x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ) ⪰ 0 , i,j ∈ Vₖ, k ∈ [p]
        for m ∈ gmpcon.make_spar_loc_con(M,t,Y)
            @constraint(model, m in PSDCone())
        end    
    # Lᵧₖ([x]ₜ₋₁ᵥₖ[x]ₜ₋₁ᵥₖᵀ ⊗ (A-xxᵀ)|_vₖ) ⪰ 0 , k ∈ [p]
        if G_con
            println("----------------G-constraints are active")
            for m ∈ gmpcon.make_spar_G_con(M,t,Y)
                @constraint(model, m in PSDCone())
            end    
        end  
    # ∑ Lᵧₖ(1)
        @objective(model, Min, gmpcon.make_obj(mc,Y))
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

end
