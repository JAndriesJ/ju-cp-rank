proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"cp_rank.jl")
include(proj_dir*"batch_run.jl")
include(proj_dir*"extract_atoms.jl")

using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .constraints ; const con = constraints
# using .cp_model
using .cp_rank
using .batch_run 
using .extract_atoms

assets_dir = dirname(dirname(proj_dir))*"\\assets\\"
# Example data
datadir = assets_dir*"data\\"
# cp_matrices.generate_random_cp_mats([5:9...], [0.25:0.075:0.95...], datadir*"\\rand\\")
# cp_matrices.gen_asc_nzd_cp_mats(9,12,datadir*"asc_nzd\\")
# cp_matrices.gen_bipart_supp_mats(5:12, datadir*"bipart\\")
# cp_matrices.generate_lit_cp_mats(datadir*"lit\\")

r_cp_mats = cp_matrices.load_mats(datadir*"rand\\")
lit_cp_mats = cp_matrices.load_mats(datadir*"lit\\")
asc_nzd_cp_mats = cp_matrices.load_mats(datadir*"asc_nzd\\")
bipart_cp_mats = cp_matrices.load_mats(datadir*"bipart\\")

t = 3
save_dir = assets_dir*"bipart\\t$(t)\\"
batch_run.batch_comp_and_save(bipart_cp_mats, t, save_dir, ["wsp","sp"]) # "id",
batch_run.make_summary(save_dir)
batch_run.clean_summary(save_dir)
batch_run.clean_clean_summary(save_dir, datadir*"bipart\\mat_data.csv")
# batch_run.get_mat_data(datadir*"bipart\\")
batch_run.plot_bounds(save_dir,t)
batch_run.plot_times(save_dir,t)


#----------------------------------------------------------
prefereti = ["ex43_n9_zd0.611_r9_ucpr23_"
             "ex51_n8_zd0.54_r8_ucpr21"
             "ex52_n9_zd0.75_r9_ucpr17"
             "ex53_n9_zd0.53_r9_ucpr25"
             "ex54_n9_zd0.64_r9_ucpr21"]

A = r_cp_mats[prefereti[5]]

A = lit_cp_mats["ex20_n8_zd0.036_r8_ucpr18_14BSUmex3_"]
A = lit_cp_mats["ex21_n11_zd0.2_r11_ucpr32_14BSUex5_"]
A = asc_nzd_cp_mats[[keys(asc_nzd_cp_mats)...][19]]

using LinearAlgebra
M, fact = cp_matrices.gen_bipart_supp_mats(9) 

A = M + diagm(rand(9))

t = 2
ξₜᶜᵖⁱᵈA  = cp_rank.get_ξₜᶜᵖ(A, t, "id") 
ξₜᶜᵖˢᵖA  = cp_rank.get_ξₜᶜᵖ(A, t, "sp")  
ξₜᶜᵖʷˢᵖA = cp_rank.get_ξₜᶜᵖ(A, t, "wsp")  

#---------------------------------------------------------
# bibartite support cp-mats




#---------------------------------------------------------
# A = r_cp_mats["ex02_n5_zd0.50_r5_ucpr10_"]
# n = size(A)[1]

# mom_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\t3\\rand\\moments\\"
# mom_path = mom_dir*"Mom_ex02_n5_zd0.50_r5_ucpr10"
# mom_path_id = mom_path*"__id_t3.csv"
# mom_path_sp = mom_path*"__sp_t3.csv"
# mom_path_wsp = mom_path*"__wsp_t3.csv"

# df_id = batch_run.load_moments(mom_path_id)
# df_sp = batch_run.load_moments(mom_path_sp)
# df_wsp = batch_run.load_moments(mom_path_wsp) 

# n_id, t_id, mom_vals_id = proc_mom(df_id)
# ext_atoms_id = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)
# cents_id, weights_id = ext_centers_weights(ext_atoms_id)
# A_ext_id = recon_mat(cents_id, weights_id)

# n_sp, t_sp, mom_vals_sp = proc_mom(df_sp)
# ext_atoms_sp =  extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-4, true)
# cents_sp, weights_sp = ext_centers_weights(ext_atoms_sp ,A)
# A_ext_sp = recon_mat(cents_sp, weights_sp)

# n_wsp, t_wsp, mom_vals_wsp = proc_mom(df_wsp)
# ext_atoms_wsp =  extract_atoms.get_atoms(n_wsp, t_wsp, mom_vals_wsp, 1e-4, true)
# cents_wsp, weights_wsp = ext_centers_weights(ext_atoms_wsp ,A)
# A_ext_wsp = recon_mat(cents_wsp, weights_wsp)

# using LinearAlgebra
# A - weights_id[1]*cents_id[1]*cents_id[1]'
# LinearAlgebra.eigvals(A - weights_id[1]*cents_id[1]*cents_id[1]')

# A
# A_ext_id 
# A_ext_sp 

# sum((A - A_ext_id).^2)
# sum((A - A_ext_sp).^2)

#----------------------------------------------------------
k = "ex02_n5_zd0.50_r5_ucpr10_"
A = r_cp_mats[k]
n = size(A)[1]

mom_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\t3\\rand\\moments\\"
mom_path = mom_dir*"Mom_"*k
mom_path_id = mom_path*"_id_t3.csv"


df_id = batch_run.load_moments(mom_path_id)

n_id, t_id, mom_vals_id = extract_atoms.proc_mom(df_id)
ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)
cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
A_ext_id                = extract_atoms.recon_mat(cents_id, weights_id)

mom_expos = moments.make_mon_expo(n,2*t_id)


nar = hcat([map(β -> sum(c.^β), mom_expos) for c in cents_id]...)
nar'*nar

nar[:,3]

c = hcat(cents_id...)*weights_id
map(β -> sum(c.^β), mom_expos)


# ####
# k = "ex02_n5_zd0.50_r5_ucpr10_"
# A = r_cp_mats[k]
# n = size(A)[1]

# mom_path = mom_dir*"Mom_"*k
mom_path_sp = mom_path*"_sp_t3.csv"
df_sp                   = batch_run.load_moments(mom_path_sp)
n_sp, t_sp, mom_vals_sp = extract_atoms.proc_mom(df_sp)
ext_atoms_sp            = extract_atoms.get_atoms(n_sp, t_sp, mom_vals_sp, 1e-5, true)
cents_sp, weights_sp    = extract_atoms.ext_centers_weights(ext_atoms_sp ,A)
A_ext_sp                = extract_atoms.recon_mat(cents_sp, weights_sp)
#####
A
A_ext_id
A_ext_sp
sqrt(sum((A - A_ext_id).^2))
sqrt(sum((A - A_ext_sp).^2))

#----------------------------------------------------------
comp_save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\"
flav = "rand" ; t = 4
comp = comp_save_dir*"\\t$(t)\\$flav\\"
dat  = datadir*"$flav\\mat_data.csv"

batch_run.make_summary(comp)
batch_run.clean_summary(comp)
batch_run.clean_clean_summary(comp, dat)

batch_run.plot_bounds(comp,t)
batch_run.plot_times(comp,t)

#----------------------------------------------------------

#Compute and save a batch----------------------------------------------------------
t = 2
save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\" # 
mkdir(save_dir)
mkdir(save_dir*"t$(t)\\")
for ex in ["rand","lit"]
    ex_save_dir = save_dir*"t$(t)\\$ex\\"
    mkdir(ex_save_dir)
    mats = cp_matrices.load_mats(datadir*"$ex\\")
    batch_run.batch_comp_and_save(mats, t,true, ex_save_dir)
    batch_run.make_summary(ex_save_dir)
end
#post process----------------------------------------------------------
t = 2
ex_save_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\"
batch_run.batch_post_proc(ex_save_dir,t)

# Flatness mass check----------------------------------------------------------
t = 4
flav = ["rand","lit"][1]
base_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\"
moment_dir = base_dir*"2022-05-10T08\\t$t\\$flav\\moments\\"
data_dir = base_dir*"data\\$flav\\"

moment_files = readdir(moment_dir)
# moment_file = moment_files[1]
# mom_path = moment_dir*moment_file

# df             = batch_run.load_moments(mom_path)
# n, t, mom_vals = extract_atoms.proc_mom(df)
# r = get_mom_mat_ranks(n, t, mom_vals)

#### ----------------------------------------------------------------------------
# moment_file = moment_files[2]
# mom_path = moment_dir .*moment_files[1]

plonk = hcat([[mf,[nar(moment_dir*mf)]] for mf in moment_files]...)

# using DataFrames, CSV
df = DataFrame( example = plonk[1,:],
                ranks   = plonk[2,:]) 
CSV.write(base_dir*"$(t)_$(flav).csv", df)


function nar(mom_path)
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
        push!(r, lowrankchchek(mom_mat, ranktol))
    end
    return r
end

data_path = data_dir*moment_files[1][5:end-10]*".csv"
A = cp_matrices.load_mats(data_path)

r = lowrankchchek(M, ranktol)
a, b, U =lowrankchol(M,1e-4)


function lowrankchol(M::AbstractMatrix, ranktol)
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
    nM, cM, (F.U[:, 1:r] * Diagonal(sqrt.(F.S[1:r])))'
end

function lowrankchchek(M::AbstractMatrix, ranktol)
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


using DataFrames, CSV
summary_csv = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\t2\\rand\\clean_clean_summary.csv"
df = CSV.read(ex_save_dir, DataFrame,delim="&")


batch_run.clean_clean_summary("C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\2022-05-10T08\\t2\\rand\\", datadir*"rand\\mat_data.csv")


# get_z_dense(M) = sum(M .== 0) / ((size(M)[1])*(size(M)[1]-1))
# vcat([ [k get_z_dense(r_cp_mats[k])] for k in keys(r_cp_mats)])

# function significant_sings(M) 
#     s = LinearAlgebra.svdvals(M)
#     return  sum(cumsum(s) .<  0.99*sum(s))
# end

# function save_load_moments(ext_mom, save_path)  #### Put this in cp_rank
#     df = DataFrame(mom_vals = ext_mom) 
#     CSV.write(save_path, df)
# end
# function save_load_moments(load_path)
#     return CSV.read(load_path, DataFrame).mom_vals
# end

#----------------------------------------------------------
# recon_mat(c::Vector{Vector{Float64}}, w) = sum([w[i]*c[i]*c[i]' for i in 1:length(c)])
# recon_mat(c::Vector{Vector{Vector{Real}}}, w) = sum([ sum([w[j][i]*c[j][i]*c[j][i]' for i in 1:length(c[j])]) for j ∈ 1:length(c)])

# function ext_centers_weights(extract ,A)
#     n = size(A)[1]
#     mc = moments.get_maximal_cliques(A) 
#     cents_sp = [[[i ∈ mc[j] ? popfirst!(a.center) : 0  for i in 1:n] for a in extract[j].atoms] for j ∈ 1:length(mc)]
#     weights_sp = [[a.weight for a in extract[j].atoms] for j ∈ 1:length(mc)]
#     return cents_sp, weights_sp
# end
# function ext_centers_weights(extract)
#     cents = [a.center for a in extract.atoms]
#     weights = [a.weight for a in extract.atoms]
#     return cents, weights
# end

# function proc_mom(df)
#     if contains(df.mom_inds[1],"]]")
#         return proc_sparse_mom(df)
#     else
#         return proc_dense_mom(df)
#     end
# end
# function proc_dense_mom(df)
#     mom_vals  = df.mom_vals 
#     mom_inds = [eval(Meta.parse(ind)) for ind in df.mom_inds]
#     n = length(mom_inds[1])
#     t = div(maximum(maximum([m for m in mom_inds])),2) 
#     return  n, t, mom_vals
# end
# function proc_sparse_mom(df)
#     mom_vals  = df.mom_vals 
#     mom_inds = [eval(Meta.parse(ind)) for ind in df.mom_inds]
#     nc = maximum([m[1] for m in mom_inds])
#     n_s = [length([m[2] for m in mom_inds if m[1] == c][1]) for c in 1:nc]
#     t = div(maximum([m[2][1] for m in mom_inds if m[1] == 1]),2) 
#     val_s = [[ mom_vals[k] for k in 1:length(mom_vals) if mom_inds[k][1] == c] for c in 1:nc]
#     return n_s, t, val_s
# end

# load_moments(load_path) = CSV.read(load_path, DataFrame)
