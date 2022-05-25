module matrix_IO

using CSV, DataFrames

### save the moment matrix ------
# save_mom_mat(ξₜᶜᵖ, n, t, save_path) = save_mat(cp_model.rec_mom_mat(n,t,ξₜᶜᵖ), save_path)
function save_moments(ξₜᶜᵖ,n,t,save_path)
    df = DataFrame(moments_expo = ξₜᶜᵖ.obj_dict[:Lx].axes[1], 
                   moments_val = cp_model.rec_mom_mat(n,2*t,ξₜᶜᵖ))
                   CSV.write(save_path, df)  
end
load_moments(save_path) = make_mom_dic(CSV.read(save_path, DataFrame))
make_mom_dic(df) = Dict(zip(clean_expo.(df.moments_expo), df.moments_val))
clean_expo(expo) = parse.(Int,map(s -> strip(s,[',','[',']']),split(expo)))

### Save the matrix ------
save_mat(M, save_path) = CSV.write(save_path, DataFrame(M, :auto))
load_mat(load_path) = Matrix(CSV.read(load_path, DataFrame))



end
