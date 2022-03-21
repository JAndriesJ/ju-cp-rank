module cp_model
using LinearAlgebra ; const la = LinearAlgebra
using JuMP
using MosekTools # The solver that we use.


proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")

using .moments ; const mom = moments
using .constraints ; const con = constraints

export modelξₜᶜᵖ,
       computeξₜᶜᵖ,
       rec_mom_mat

function modelξₜᶜᵖ(A,t,conlist)
    n = size(A)[1]
    model = Model()
    @variable(model, Lx[mom.make_mon_expo(n,2*t)] ) # Define variables in the moment matrix.
## Build the moment matrix and constrain it to be PSD.
    mom_matₜ_expo = mom.make_mon_expo(n,(t,t),A)
    mom_matₜ      = con.α_to_Lxᵅ(Lx, mom_matₜ_expo)
    @constraint(model, la.Symmetric(mom_matₜ) in PSDCone())
# Second order Moment constraints
    mom_mat₌₁ = mom.make_mon_expo(n,(1,1),isle=false)
    Lx_mom_mat₌₁ = con.α_to_Lxᵅ(Lx,mom_mat₌₁)
    for  i in 1:n, j in i:n
        fix(Lx_mom_mat₌₁[i,j], A[i,j])
    end
# Localizing g constraint
    for c in con.make_loc_con(A,t,Lx)
        @constraint(model, c in PSDCone())
    end
# Dagger constraints
    if occursin("Dag",conlist)
        dag_con  = con.make_dag_con(A,t,Lx)
        for key in keys(dag_con)
            @constraint(model, dag_con[key] .>= zeros(size(dag_con[key])))
        end
        println("----------------Dagger constraints are active")
    end
# Localizing XX constraint
    if occursin("XX",conlist)
        println("----------------XX constraints are active")
        for c in con.make_xx_con(A,t,Lx)
            @constraint(model, c in PSDCone() )
        end
    end
# G Constraints weak and strong
    if occursin("sG",conlist)
        println("----------------G-constraints are active")
        G_con = con.make_G_con(A,t,Lx)
        @constraint(model, la.Symmetric(G_con) in PSDCone())
    end
# Ideal constraints
    if occursin("id",conlist)
        println("----------------G-constraints are active")
        for c in con.make_ideal_con(A,t,Lx)
            @constraint(model, c .== 0)
        end
    end


##  Set objective
    @objective(model, Min, mom_matₜ[1,1])
    return model
end


"""This is where the ξₜᶜᵖ is calculated for matrix A """
function computeξₜᶜᵖ(model)
    set_optimizer(model, Mosek.Optimizer)
    optimize!(model) 
    # output results
    println("Primal: ", primal_status(model))
    println("Dual: ", dual_status(model))
    println("Objective: ", objective_value(model))
    return model
end



#Post processing--------------------------------------------------------------

function rec_mom_mat(n::Int64,t::Tuple{Int64, Int64},ξₜᶜᵖ)
    Lx      = ξₜᶜᵖ.obj_dict[:Lx]
    MB_exp  = mom.make_mon_expo(n,t)
    MB      = con.α_to_Lxᵅ(Lx, MB_exp)
    return JuMP.value.(MB)
end

function rec_mom_mat(n::Int64,t::Int64,ξₜᶜᵖ)
    Lx      = ξₜᶜᵖ.obj_dict[:Lx]
    MB_exp  = mom.make_mon_expo(n,t)
    MB      = con.α_to_Lxᵅ(Lx, MB_exp)
    return JuMP.value.(MB)
end

function rec_mom_mat(A,ξₜᶜᵖ)
    Lx      = ξₜᶜᵖ.obj_dict[:Lx]
    MB      = con.α_to_Lxᵅ(Lx,A)
    return JuMP.value.(MB)
end

function rec_mom_mat(n::Int64,t::Int64, ξₜᶜᵖ, mon_cliques)
    mon = mom.make_mon_expo(n,t)
    cliq_ordering = mom.order_monomials(mon, mon_cliques)
    cliq_ordered_mon = mom.make_mon_expo(mon[cliq_ordering])
    return rec_mom_mat(cliq_ordered_mon, ξₜᶜᵖ)
end




function check_faltness(n::Int64,t1::Int64,t2::Int64, ξₜ₁ᶜᵖ, ξₜ₂ᶜᵖ)
    mom_mat_1 = rec_mom_mat(n,t1,ξₜ₁ᶜᵖ)
    mom_mat_2 = rec_mom_mat(n,t2,ξₜ₂ᶜᵖ)
    return la.rank(mom_mat_1) == la.rank(mom_mat_2)
end


end


# if occursin("wG",conlist)
#     println("----------------Weak G-constraints are active")
#     weakG_con = con.make_weakG_con(A,t,Lx)
#     for key in keys(weakG_con)
#         weakG_con_key = weakG_con[key]
#         m = size(weakG_con_key)[1]
#        @constraint(model, la.Symmetric(weakG_con_key) in PSDCone())
#     end
# end