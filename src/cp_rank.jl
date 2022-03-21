module cp_rank

using CSV, DataFrames

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
# include(proj_dir*"moments.jl")
# include(proj_dir*"constraints.jl")
include(proj_dir*"cp_model.jl")

using .cp_matrices ; const cpm = cp_matrices
# using .moments ; const mom = moments
# using .constraints ; const con = constraints
using .cp_model  


export get_ξₜᶜᵖ,
       save_mat,
       load_mat,
       save_mom_mat,
       run_and_save_get_ξₜᶜᵖ,
       get_Bomze_cp_mats,
       get_random_cp_mat,
       get_random_sparse_cp_mat,
       gen_random_band_mat
           
get_Bomze_cp_mats() = cpm.get_Bomze_cp_mats()
get_random_cp_mat(n,r) = cpm.get_random_cp_mats(n,r)
get_random_sparse_cp_mat(M::Matrix{Int64}) = cpm.get_random_sparse_cp_mats(M)
get_random_sparse_cp_mat(n::Int,p::Float64=0.5) = cpm.get_random_sparse_cp_mats(cpm.gen_random_sparcity_mat(n,p))
gen_random_band_mat(n,k) = cpm.gen_random_band_mat(n,k)


get_ξₜᶜᵖ(M,t,conlist) = cp_model.computeξₜᶜᵖ(cp_model.modelξₜᶜᵖ(M,t,conlist))



### Save the computations ------
function run_and_save_get_ξₜᶜᵖ(rsM_n ,t ,cons, save_path)
    original_stdout = stdout;
    (rd, wr) = redirect_stdout();
    ξₜᶜᵖrsM_n_t_cons = get_ξₜᶜᵖ(rsM_n ,t ,cons)
    s = []
    for i in 1:100
        rl = readline(rd)
        push!(s,rl)
        if contains(rl,"Objective:")
            break
        end
    end
    redirect_stdout(original_stdout);

    ###---
    touch(save_path)
    open(save_path, "w") do io
        for line in s
            write(io, line*"\n")
        end
    end;
    return ξₜᶜᵖrsM_n_t_cons
end

### save the moment matrix ------
save_mom_mat(ξₜᶜᵖ, n, t, save_path) = save_mat(cp_model.rec_mom_mat(n,t,ξₜᶜᵖ), save_path)

### Save the matrix ------
save_mat(M, save_path) = CSV.write(save_path, DataFrame(M, :auto))
load_mat(load_path) = Matrix(CSV.read(load_path,DataFrame))


function run_tests()
    cpm.run_tests()
    mom.run_tests()
    con.run_tests()
    cp_model.run_tests()  
    compute.run_tests() 
end

end
