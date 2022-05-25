proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"cp_rank.jl")

using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .constraints ; const con = constraints
using .cp_model
using .cp_rank

include(proj_dir*"batch_run.jl")
using .batch_run 


# Load matrices -------------------------------------------------------------------
## From literature
lit_cp_mats = cp_matrices.get_lit_cp_mats()# load a dictionary of cp matrices from literature.
M = lit_cp_mats[6]
cp_rank.save_mat(M,load_path)
## Randomly generated
random_mats_path = dirname(dirname(proj_dir))*"\\assets\\data\\randomly_generated_cp-matrices\\mats\\"
matrix_names = readdir(random_mats_path)
M_path = random_mats_path*matrix_names[22]
M = cp_rank.load_mat(load_path)
## Non-examples


# Visualizations
cp_rank.show_mat_support(M) 
cp_rank.show_support_graph(M)
ze = mom.get_zero_entries(M)
nze = mom.get_nonzero_entries(M)
t = 2
cp_rank.show_mat_support(mom.get_mom_mat_supp(t,M))
cp_rank.show_mat_support(mom.get_ten_con_supp(t,M))
## Compute ξₜᶜᵖ(M) ---------------------------------------------
ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M, t, "sG")
ξₜᶜᵖid = cp_rank.get_ξₜᶜᵖ(M, t, "sGid")

# Compute ξₜᶜᵖ(M) and save output
# save_dir = "C:\\Users\\jandr\\Dropbox\\Vaults\\Pandemonium\\LABOURS\\4-Toulouse\\Numerics\\"*string(now())[1:10]*"\\"
# mkdir(save_dir)
cp_rank.run_and_save_get_ξₜᶜᵖ(M,t,conlist,save_dir*"ξₜᶜᵖrsM_n$(n)_t$(t)_$conlist.txt")


## --------------------------------------------
## --------------------------------------------
## --------------------------------------------
## --------------------------------------------

## Sandbox part
assets_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\data\\randomly_generated_cp-matrices\\"
save_dir = assets_dir*"mats\\"
# batch_run.gen_batch_sparse_cp_mats(5:9,1:10,save_dir)
load_dir = assets_dir*"mats\\"
Mats_dict = batch_run.load_batch_sparse_cp_mats(load_dir)
mat_name = "M_n7_ex02_r7_rp16"
M = Mats_dict[mat_name]


batch_run.get_batch_mat_support(load_dir,3)
# band matrices
for k = 1:8
    M = cp_rank.gen_random_band_mat(8,k)
    mat_supp = cp_rank.show_mat_support(M);
    t = 2 
    mom_mat_supp = cp_rank.show_mat_support(moments.get_mom_mat_supp(t,M));
    ten_con_supp = cp_rank.show_mat_support(moments.get_ten_con_supp(t,M)); 
    using FileIO

    save("supp_k$(k).png", mat_supp)
    save("mom_supp_k$(k)_t$t.png", mom_mat_supp)
    save("ten_supp_k$(k)_t$t.png", ten_con_supp)
end

# using CSV, DataFrames
# using Missings
# clean_csv_path = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022.03.16\\clean_summary.csv"
# clean_csv = CSV.read(clean_csv_path,DataFrame)
# mat_sizes = Matrix(clean_csv[:,[ :sm, :sms, :smg, :smgs]])

# using Plots.PlotMeasures
# scatter([1:27], mat_sizes[:,[1,2]],xticks=(1:27,readdir(load_dir)[1:27]),xrotation=90, bottom_margin = 130px)
# scatter([1:27], mat_sizes[:,[3,4]],xticks=(1:27,readdir(load_dir)[1:27]),xrotation=90, bottom_margin = 130px)


# comps = collect(Missings.replace(Matrix(clean_csv[:,[ :sG_obj_v, :sGid_obj_v, :sG_time_s, :sGid_time_s]]), 0.0))
# comps_spar = map(a -> parse(Float64,a), replace(comps,"-" => "0.0")[:,[2,4]])
# comps_dens = comps[:,[1,3]]

# scatter([1:27], [comps_spar[:,1], comps_dens[:,1]], xticks=(1:27,readdir(load_dir)[1:27]),xrotation=90, bottom_margin = 130px,legend = false) # objvals
# scatter([1:27], [comps_spar[:,2], comps_dens[:,2]], xticks=(1:27,readdir(load_dir)[1:27]),xrotation=90, bottom_margin = 130px,legend = false) # times 


conlist = "sGid" #a string that if it contains sG it impliments the G⊗L ⪰ 0 constriants, If it contains id it impliments the ideal constraints.
t = 2
ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,conlist)

comp_dir = assets_dir*"comp\\"
# batch_run.batch_comp_and_save(Mats_dict, 3, comp_dir)
# Compute rank of moment matrix
# Make larg scale vizualizations.
assets_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\16.03\\"
batch_run.make_summary(assets_dir)
batch_run.clean_summary(assets_dir)


for fmat_name in readdir(load_dir)
    mat_name = fmat_name[1:end-4]
    M = Mats_dict[mat_name]
    G_plot = cp_rank.show_support_graph(M)
    save(load_dir*"supp_G_"*mat_name*".png",G_plot)
    # max_cli = moments.get_maximal_cliques(M)
    # cp_rank.save_mat([max_cli] ,load_dir*"max_cli_"*mat_name*".csv")
end

















