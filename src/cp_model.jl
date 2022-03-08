module cp_model
using LinearAlgebra
const la = LinearAlgebra
using JuMP


proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
using .moments ; const mom = moments
using .constraints ; const con = constraints

export modelξₜᶜᵖ,
       rec_mom_mat

function modelξₜᶜᵖ(A,t,conlist)
    n = size(A)[1]
    model = Model()
    @variable(model, Lx[mom.make_mon_expo(n,2*t)] ) # Define variables in the moment matrix.
## Build the moment matrix and constrain it to be PSD.
    mom_matₜ_expo = mom.make_mon_expo(n,(t,t))
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

function rec_mom_mat(n::Int64,t::Int64,ξₜᶜᵖ)
    Lx      = ξₜᶜᵖ.obj_dict[:Lx]
    MB_exp  = mom.make_mon_expo(n,(t,t))
    MB      = con.α_to_Lxᵅ(Lx, MB_exp)
    return JuMP.value.(MB)
end
function rec_mom_mat(A::Matrix{Vector{Int64}},ξₜᶜᵖ)
    Lx      = ξₜᶜᵖ.obj_dict[:Lx]
    MB      = con.α_to_Lxᵅ(Lx,A)
    return JuMP.value.(MB)
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