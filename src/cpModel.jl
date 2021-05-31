module cpModel
using LinearAlgebra
const la = LinearAlgebra
using JuMP


proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"Utility.jl")
include(proj_dir*"Moments.jl")
include(proj_dir*"Constraints.jl")
using .Utility
using .Moments
using .Constraints

export Modelξₜᶜᵖ

function Modelξₜᶜᵖ(A,t,conlist)
    n = size(A)[1]
    model = Model()
    list_of_keys = make_mom_expo_keys(n, t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] )
## Build the moment matrix and constrain it to be PSD.
    mom_matₜ_expo = Moments.make_mon_expo_mat(n,(t,t),true)
    mom_matₜ      = Utility.index_to_var(Lx, mom_matₜ_expo)
    @constraint(model, la.Symmetric(mom_matₜ) in PSDCone())
# Second order Moment constraints
    mom_mat₌₁ = Moments.make_mon_expo_mat(n,(1,1),false)
    Lx_mom_mat₌₁ = Utility.index_to_var(Lx,mom_mat₌₁)
    for  i in 1:n, j in i:n
        fix(Lx_mom_mat₌₁[i,j], A[i,j])
    end
# Localizing g constraint
    loc_con = Constraints.make_loc_con(A,t,Lx)
    Z_mat = zeros(size(loc_con[(1,1)]))
    for key in keys(loc_con)
        @SDconstraint(model, loc_con[key] >= Z_mat)
    end
# Dagger constraints
    if occursin("Dag",conlist)
        dag_con  = Constraints.make_dag_con(A,t,Lx)
        for key in keys(dag_con)
            @constraint(model, dag_con[key] .>= zeros(size(dag_con[key])))
        end
        println("----------------Dagger constraints are active")
    end

# Localizing XX constraint
    if occursin("XX",conlist)
        println("----------------XX constraints are active")
        xx_con = Constraints.make_xx_con(A,t,Lx)
        for key in keys(xx_con)
            @SDconstraint(model, xx_con[key] >= zeros(size(xx_con[key])))
        end
    end
# G Constraints weak and strong
    if occursin("wG",conlist)
        println("----------------Weak G-constraints are active")
        weakG_con = Constraints.make_weakG_con(A,t,Lx)
        for key in keys(weakG_con)
            weakG_con_key = weakG_con[key]
            m = size(weakG_con_key)[1]
           @constraint(model, la.Symmetric(weakG_con_key) in PSDCone())
        end
    end
    if occursin("sG",conlist)
        println("----------------G-constraints are active")
        G_con                 = Constraints.make_G_con(A,t,Lx)
        @constraint(model, la.Symmetric(G_con) in PSDCone())
    end

##  Set objective
    @objective(model, Min, Lx[zeros(n)])
    return model
end



end
