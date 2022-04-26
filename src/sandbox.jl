proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
include(proj_dir*"GMP_constraints.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"cp_rank.jl")

using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .constraints ; const con = constraints
using .GMP_constraints ; const gmpcon = GMP_constraints
using .cp_model
using .cp_rank

include(proj_dir*"batch_run.jl")
using .batch_run 
#----------------------------------------------------------
datadir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\data\\randomly_generated_cp-matrices\\mats\\"
r_cp_mats = cp_matrices.get_random_cp_mats(datadir)
save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.04.26\\t2\\"
batch_run.batch_comp_and_save(r_cp_mats, 2,true, save_dir)
batch_run.make_summary(save_dir)
batch_run.clean_summary(save_dir)
#----------------------------------------------------------
lit_cp_mats = cp_matrices.get_lit_cp_mats()
save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.04.25\\lit\\t3\\"
batch_run.batch_comp_and_save(lit_cp_mats, 3,true, save_dir)

batch_run.make_summary(save_dir)
batch_run.clean_summary(save_dir)

#----------------------------------------------------------
# M = r_cp_mats["M_n5_ex02_r5_rp9"]
# t = 2
# idea = cp_rank.get_ξₜᶜᵖ(M, t, "ideal", true) 
# gmpd = cp_rank.get_ξₜᶜᵖ(M, t, "gmpd", true)  # slow
# gmps = cp_rank.get_ξₜᶜᵖ(M, t, "gmps", true)  #wrong
#----------------------------------------------------------






#### result Visualizations
# using Plots
# using Plots.PlotMeasures
# using CSV, DataFrames
# using FileIO
# sum_path = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.04.11\\t3\\matrix localizing constraints\\clean_summary.csv"
# sum_path = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.04.11\\t3\\no matrix localizing constraints\\clean_summary.csv"
# df = CSV.read(sum_path,DataFrame,delim =",")
# n = size(df)[1]

# objs = scatter([1:n], Matrix(df[:,[:ideal_obj_v, :gmpd_obj_v, :gmps_obj_v]]),
#                     xticks=(1:n,df.name),
#                     xlabel="example",
#                     ylabel="bound",
#                     xrotation=90,
#                     bottom_margin = 130px,
#                     label = ["ξₜᶜᵖ⁻ⁱᵈ" "ξₜᶜᵖ⁻ᴳᴹᴾ-dense" "ξₜᶜᵖ⁻ᴳᴹᴾ-sparse"],
#                     legend=:bottomright,
#                     fillalpha = 0.1,
#                     markershape = [:circle :star5 :diamond])                  
# img_save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.04.11\\"
# save(img_save_dir*"t3_no_matrix_loc.png", objs)

# cp_matrices.run_tests()
# moments.run_tests()



