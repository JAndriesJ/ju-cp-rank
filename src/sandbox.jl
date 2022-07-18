#-------------------------------------Diagnostics--------------------------------------------
test_dir = dirname(dirname(@__FILE__))*"\\test\\"
include(test_dir*"matrix_IO_test.jl")
include(test_dir*"moments_test.jl")
include(test_dir*"cp_matrices_test.jl")
include(test_dir*"cp_model_test.jl")

#-------------------------------------Initialization--------------------------------------------
src_dir = dirname(@__FILE__)*"\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_matrices.jl")
include(src_dir*"cp_model.jl")
include(src_dir*"batch_run.jl")

# include(proj_dir*"extract_atoms.jl")
#-------------------------------------Data--------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
exa = ["rand" "lit" "nlit"][2]*"\\"
dataexadir = data_dir*exa
#-------------------------------------Data--------------------------------------------
# cp_matrices.generate_random_cp_mats(dataexadir)
# cp_matrices.generate_lit_cp_mats(dataexadir)
# cp_matrices.generate_lit_non_cp_mats(dataexadir)

#-------------------------------------Meta Data--------------------------------------------
# matrix_IO.get_mat_data(dataexadir)

#-------------------------------------Batch computations--------------------------------------------#
mats = matrix_IO.load_mats(dataexadir)
t = 2
results_dir = assets_dir*"results\\$exa"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_subdir = results_dir*"t$(t)\\"
!isdir(results_subdir) ? mkdir(results_subdir) : 0

# batch_run.batch_comp_and_save(mats, t, "G".*["wsp", "sp", "id"], results_subdir) # "",  .*"Gddagxx", 
batch_run.batch_comp_and_save(mats, t, "Gdag".*["wsp", "sp", "id"], results_subdir)
# batch_run.batch_comp_and_save(mats, t, "Gddagxx".*["wsp", "sp", "id"], results_subdir)

batch_run.make_summary(results_subdir)
batch_run.clean_summary(results_subdir, dataexadir)
batch_run.plot_times(results_subdir, t)
#------------------------------------Post Processing-------------------------------------------------
batch_run.make_summary(results_subdir)
batch_run.clean_summary(results_subdir)
batch_run.clean_clean_summary(results_subdir, datadir*"$exa\\mat_data.csv")

batch_run.plot_bounds(results_subdir, t)
batch_run.plot_times(results_subdir, t)
#-------------------------------------------------------------------------
moment_dir = assets_dir*"results\\$exa\\t$(t)\\moments\\"
batch_run.get_all_mom_mat_ranks(moment_dir)

#----------------------------Running specific matrices---------------------------------------------
M = mats[[keys(mats)...][20]]
t = 1

ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"");
ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"G");
ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gdag");
ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gddag");
ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gddagxx");

ξₜᶜᵖʷˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"");
ξₜᶜᵖʷˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"G");
ξₜᶜᵖʷˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gdag");
ξₜᶜᵖʷˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gddag");
ξₜᶜᵖʷˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gddagxx");

ξₜᶜᵖˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"");
ξₜᶜᵖˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"G");
ξₜᶜᵖˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gdag");
ξₜᶜᵖˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gddag");
ξₜᶜᵖˢᵖ , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gddagxx");

#------------------------------------Graphs-------------------------------------
using Graphs
using GraphPlot, Cairo, Fontconfig
G = Graph(M  .> 0)
gplot(G,nodelabel=1:size(M)[1])
#-------------------------------------------------------------------------

#----------------------------------------------------------
# mom_dir = assets_dir*"asc_nzd_6\\t3\\moments\\" 
# mom_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\asc\\old\\t3\\moments\\"
mom_dir = assets_dir*"220706\\2022-05-10T08\\t3\\rand\\moments\\"
mom_dir = assets_dir*"results\\lit\\t3\\moments\\"
mom_dir_list = [c for c in readdir(mom_dir) if contains(c,"_id")]

mom_path = mom_dir*mom_dir_list[3] # 2 works
df_id = matrix_IO.load_moments(mom_path)
n_id, t_id, mom_vals_id = extract_atoms.proc_mom(df_id)
ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)
cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
A_ext_id                = extract_atoms.recon_mat(cents_id, weights_id)
length(weights_id)

# A =  matrix_IO.load_mat(assets_dir*"data\\Xold\\rand\\ex02_n5_zd0.50_r5_ucpr10_.csv") 
A =  matrix_IO.load_mat(assets_dir*"data\\lit\\ex03_n5_nzd0.5_r5_ucpr5_98Bex4.3_.csv") 
mom_path_sp = mom_dir*[c for c in readdir(mom_dir) if contains(c,"_sp")][3]
df_sp                   = batch_run.load_moments(mom_path_sp)
n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(df_sp)
ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4, true)

cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp ,A)
A_ext_sp                = extract_atoms.recon_mat(cents_sp, weights_sp)


# Flatness mass check----------------------------------------------------------
using LinearAlgebra
using DataFrames, CSV

moment_dir = assets_dir*"results\\$exa\\t2\\moments\\"
get_all_mom_mat_ranks(moment_dir )

function get_all_mom_mat_ranks(moment_dir::String)
    moment_files = readdir(moment_dir)
    mat_df = hcat([[mf,[get_mom_mat_ranks(moment_dir*mf)]] for mf in moment_files]...)
    df = DataFrame( example = mat_df[1,:],
                    ranks   = mat_df[2,:]) 
    CSV.write(dirname(moment_dir)*"_ranks.csv", df)
end

function get_mom_mat_ranks(mom_path::String)
    df             = batch_run.load_moments(mom_path)
    n, t, mom_vals = extract_atoms.proc_mom(df)
    return get_mom_mat_ranks(n, t, mom_vals)
end
get_mom_mat_ranks(n::Vector{Int64}, t, mom_vals) = [get_mom_mat_ranks(n[j], t, mom_vals[j]) for j in 1:length(n)  ]
function get_mom_mat_ranks(n::Int64, t, mom_vals)
    r = []
    for j in 1:t
        mom_mat_thing  = extract_atoms.get_mom_mat(n, j, mom_vals)
        mom_mat        = [m for m in mom_mat_thing.Q]
        push!(r, lowrankchchek(mom_mat))
    end
    return r
end

function lowrankchchek(M::AbstractMatrix, ranktol =1e-4)
    F = svd(M)
    nM = F.S[1] # norm of M
    tol = nM * ranktol
    r = something(findfirst(σ2 -> σ2 <= tol, F.S), 0)
    if iszero(r)
        cM = ranktol
        r = length(F.S)
    else
        cM = F.S[r] / nM
        r -= 1
    end
    r
end


get_clique_rank(mats["ex02_n5_nzd0.5_r5_ucpr0_21BS-Mp302_"])
using LinearAlgebra
get_clique_rank.([mats[k] for k in keys(mats)])
function get_clique_rank(M)
    n = size(M)[1]
    mc = moments.get_maximal_cliques(M)
    E = [[i,j] for i ∈ 1:n, j ∈ 1:n if (j > i && M[i,j] != 0 )]
    CE_mat = [Set(e) ⊆ Set(c) for e in E, c in mc]
    return rank(CE_mat)
end



