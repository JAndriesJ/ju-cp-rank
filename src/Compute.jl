module compute
using JuMP # For the optimization frame work.
using MosekTools # The solver that we use.

export computeξₜᶜᵖ,
       rec_mom_mat

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

function rec_mom_mat(n::Int64,t::Int64,Lx)
    MB_exp  = make_mon_expo_mat(n,t)
    MB      = index_to_var(Lx, MB_exp)
    mom_mat = value.(MB)
    return mom_mat
end
function rec_mom_mat(A::Array{Float64,2},t::Int64,Lx)
    @assert size(A)[1] == size(A)[2]
    n = size(A)[1]
    return rec_mom_mat(n::Int64,t::Int64,Lx)
end

end  



# extract_model_stats(model_opt) = (string(primal_status(model_opt)), string(dual_status(model_opt)), objective_value(model_opt))