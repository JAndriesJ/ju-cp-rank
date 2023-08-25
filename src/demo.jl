#-------------------------------------Diagnostics--------------------------------------------
include(dirname(dirname(@__FILE__))*"\\test\\"*"runtests.jl")

#-------------------------------------Initialization--------------------------------------------
src_dir = dirname(@__FILE__)*"\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")

include(src_dir*"cp_matrices.jl")
include(src_dir*"cp_model.jl")

# ------------------------------------Completely positive matrix rank bounds---------------------------------------------------

M =    [1.0       0.707107  0.0       0.0       0.447214
        0.707107  1.0       0.408248  0.0       0.0
        0.0       0.408248  1.0       0.288675  0.0
        0.0       0.0       0.288675  1.0       0.223607
        0.447214  0.0       0.0       0.223607  1.0]
lvl = 1
cons  = join(["G","dag","ddag","xx"][[1,2,4]],"") # choose which constraitns to inclued (see ... in the publication)
hier = ["id", "sp", "wsp"][2]          # choose from the three main hierarchies (see ... in the publication)

# Dense completely positive moment hierarchy
ξₜᶜᵖⁱᵈ , moment_matrix = cp_model.get_ξₜᶜᵖ(M,lvl,hier*cons);        





#-------------------------------------Load Data--------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
mat_types = ["Euclidean_distance_matrices"
             "bipartite"
             "literature"
             "literature_alt"
             "random_completely_positive"
             "random_nonnegative"]

mat_dir = data_dir*mat_types[3]*"\\"
data_matrix_dict = matrix_IO.load_mats(mat_dir)
dmdk = data_matrix_dict.keys
dmdv = data_matrix_dict.vals

k = 1
M_name = dmdk[k]
M      = dmdv[k]

#--------------------------ATOM EXTRACTION--------------------------------
mom_dir = assets_dir*"results\\lit\\t2\\moments\\"
mom_dir_list = [c for c in readdir(mom_dir) if contains(c,"id")]
mom_path = mom_dir*mom_dir_list[1] 

df_id = matrix_IO.load_moments(mom_path)
mom_path


include(src_dir*"extract_atoms.jl")


n_id, t_id, mom_vals_id = extract_atoms.proc_mom(df_id)
ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)

cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)

A_ext_id                = extract_atoms.recon_mat(cents_id, weights_id)


# A =  matrix_IO.load_mat(assets_dir*"data\\Xold\\rand\\ex02_n5_zd0.50_r5_ucpr10_.csv") 
A =  matrix_IO.load_mat(assets_dir*"data\\lit\\ex03_n5_nzd0.5_r5_ucpr5_98Bex4.3_.csv") 
mom_path_sp = mom_dir*[c for c in readdir(mom_dir) if contains(c,"_sp")][3]

df_sp                   = batch_run.load_moments(mom_path_sp)
n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(df_sp)
ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4, true)

cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp ,A)
A_ext_sp                = extract_atoms.recon_mat(cents_sp, weights_sp)



# ------------------------------------Nonnegative---------------------------------------------------
# ------------------------------------Nonnegative---------------------------------------------------


# 







#-------------------------------------Generating Data--------------------------------------------
# cp_matrices.generate_random_cp_mats(dataexadir)
# cp_matrices.generate_lit_cp_mats(dataexadir)
# cp_matrices.generate_lit_non_cp_mats(dataexadir)
# cp_matrices.gen_bipart_supp_mats([[2,3],[3,3],[3,4],[4,4],[4,5],[5,5],[5,6],[6,6],[6,7],[7,7],[7,8],[8,8]], dataexadir)


#-------------------------------------Batch computations--------------------------------------------#
include(src_dir*"batch_run.jl")
mats = matrix_IO.load_mats(dataexadir)
t = 3
results_dir = assets_dir*"results\\$exa"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_subdir = results_dir*"t$(t)\\"
!isdir(results_subdir) ? mkdir(results_subdir) : 0

batch_run.batch_comp_and_save(mats, t, "cp G ddag xx".*["wsp", "sp", "id"], results_subdir) # "",  .*"Gddagxx", 
# batch_run.batch_comp_and_save(mats, t, "cp Gdag".*["wsp"], results_subdir)
# batch_run.batch_comp_and_save(mats, t, "cp Gdag".*["sp"], results_subdir)
# batch_run.batch_comp_and_save(mats, t, "cp Gdag".*["id"], results_subdir)
# batch_run.batch_comp_and_save(mats, t, "Gddagxx".*["wsp", "sp", "id"], results_subdir)
#------------------------------------Post Processing-------------------------------------------------
batch_run.make_summary(results_subdir)
batch_run.clean_summary(results_subdir, dataexadir)
# batch_run.plot_bounds(results_subdir, t)
batch_run.plot_times(results_subdir, t)




using CSV, DataFrames
using Plots
using Plots.PlotMeasures
df = CSV.read(results_subdir*"clean_summary.csv", DataFrame,delim =",")
sort!(df, [:name])
n = length(df.n)
times = names(df)[end-2:end]
# df[:,["n","nzd"]]
nnzd = df[:,["n","nzd"]]
poes = [[a.n,a.nzd] for a ∈ eachrow(nnzd)]
narara = unique([ [a,b] for (a,b) in eachrow(nnzd)])
p_dict = Dict(zip(narara,[1:length(narara)...]))
nar =  map(p -> p_dict[p],poes)
n = length(narara)
dara = deepcopy(Matrix(df[:,times]))

for i ∈ 1:849, j ∈ 1:3
    if ismissing( dara[i,j])
        continue
    elseif dara[i,j] > 1100
        dara[i,j] = missing
    end
end

objs = scatter(nar , dara .+ 1.0,  # , :ucpr
                                    xticks=(1:n, [  a[2]==1 ? (Int(a[1]),a[2]) : a[2]  for a ∈ narara]),#df.ex .*"_nzd=".* string.(df.nzd)),# xlabel="Non-zero density for $(np)x$(np) matrices",
                                    yticks= [1,10,100,1000],
                                    xlabel="Size and nonzero density of the matrix",
                                    ylabel="time(s)",
                                    xrotation=65,
                                    left_margin = 20px,
                                    bottom_margin = 60px,
                                    label = ["id" "sp" "wsp"], # "e-rank_cp"
                                    legend=:topleft,
                                    fillalpha = 0.1,
                                    yscale=:log10,
                                    ylims= (0.75, 9500),
                                    markershape = [:square :diamond :circle],
                                    markersizes = [4 4 3],
                                    markercolors = [:red :yellow :green],
                                    title = "Hierarchy times at t=$t",
                                    size = (1200, 600) )     #  




#-------------------------------------------------------------------------
moment_dir = assets_dir*"results\\$exa\\t$(t)\\moments\\"
batch_run.get_all_mom_mat_ranks(moment_dir)



