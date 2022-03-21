# using Random
# Random.seed!(373)

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

# using Dates

## Demo part
# Load matrices -------------------------------------------------------------------
# Bomze_mats_dict = cp_rank.get_Bomze_cp_mats() # load a dictionary of dense cp matrices from literature.
# M = Bomze_mats_dict["M7t"] # pick a matrix from the dictionary
M = cp_matrices.get_Nie_Mats()[2]
# Generate matrices
n = 8 # set a size for matrix generation
r = 5 # set upper bound for the cp-rank of the matrix
random_cp_mat = cp_rank.get_random_cp_mat(n,r)
p = 0.3 # set a probability
random_sparse_cp_mat = cp_rank.get_random_sparse_cp_mat(n,p)
k = 3 # set a band width ≤ 4
random_band_mat = cp_rank.gen_random_band_mat(n,k)
# Load/save a pregenerated matrix
load_path = "C:\\Users\\jandr\\Dropbox\\Vaults\\Pandemonium\\LABOURS\\4-Toulouse\\Numerics\\2022-03-07T14-52-40.318\\cp_mats_6\\Ex_3rsM_6.csv"
M = cp_rank.load_mat(load_path)
cp_rank.save_mat(M,load_path)
#vizualize the support of the matrices ---------------------------------------------
cp_rank.show_mat_support(random_sparse_cp_mat) 
cp_rank.show_mat_support(random_band_mat) 
## Compute ξₜᶜᵖ(M) ---------------------------------------------
r,R,M = random_sparse_cp_mat
t = 5 # level of the hierarchy
conlist = "sGid" #a string that if it contains sG it impliments the G⊗L ⪰ 0 constriants, If it contains id it impliments the ideal constraints.
ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,conlist)
# Compute ξₜᶜᵖ(M) and save output
save_dir = "C:\\Users\\jandr\\Dropbox\\Vaults\\Pandemonium\\LABOURS\\4-Toulouse\\Numerics\\"*string(now())[1:10]*"\\"
mkdir(save_dir)
cp_rank.run_and_save_get_ξₜᶜᵖ(M,t,conlist,save_dir*"ξₜᶜᵖrsM_n$(n)_t$(t)_$conlist.txt")
## --------------------------------------------

sparsity.show_mat_support(M)
sparsity.show_support_graph(M)
# ze = sparsity.get_zero_entries(M)
# nze = sparsity.get_nonzero_entries(M)

## Sandbox part
include(proj_dir*"batch_run.jl")
using .batch_run 

assets_dir = "C:\\Users\\jandr\\code_projects\\ju-cp-rank\\assets\\16.03\\"
save_dir = assets_dir*"mats\\"
# batch_run.gen_batch_sparse_cp_mats(5:9,1:10,save_dir)
load_dir = assets_dir*"mats\\"
Mats_dict = batch_run.load_batch_sparse_cp_mats(load_dir)
M = Mats_dict["M_n7_ex02_r7_rp16"]
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



using LinearAlgebra
matrix_name = "M_n6_ex09_r6_rp8"
M = Mats_dict[matrix_name]
u, s, v  = svd(M)


mom_dict = batch_run.load_moments(assets_dir*"moments\\Moments_"*matrix_name*"_t3_sGid.csv")
mom_mat_val_seq = batch_run.get_mom_mat_val_seq(6,3,mom_dict)
Kat = mom_mat_val_seq[1]

u, s, v  = svd(mom_mat_val_seq[1])
1 + sum((s[2:end] ./ s[1:end-1]) .≥ 1e-3)
rank(mom_mat_val_seq[3])

using Plots
plot(s)



plot(cumsum(s))
sum(s .> sum(s)/length(s))
cumsum(s)

# function batch_get_matrix_sizes(loaddir,t)
load_dir = assets_dir*"Mats\\"
Mats_dict = batch_run.load_batch_sparse_cp_mats(load_dir)
for mat in readdir(load_dir)
    M = Mats_dict[mat[1:end-4]]
    siz_Mₜ, siz_Mₜ_spar, siz_MₜG, siz_MₜG_spar = get_matrix_sizes(M,t)
end


function get_matrix_sizes(M,t)
    n = size(M)[2]
    siz_Mₜ = length(mom.make_mon_expo(n,t,M))
    siz_Mₜ_spar =length(mom.make_mon_expo(n,t))
    siz_MₜG =length(mom.make_mon_expo(n,t-1,M))*n
    siz_MₜG_spar =length(mom.make_mon_expo(n,t-1))*n
    return siz_Mₜ, siz_Mₜ_spar, siz_MₜG, siz_MₜG_spar
end



# mc = sort(sparsity.get_maximal_cliques(M))
mon_cliques = sparsity.get_monomial_cliques(n,t,M)
spar_sup = mom.get_monomial_suport(n,t,M)

for t = 0:5, n = 4:5
    M = cp_rank.get_random_sparse_cp_mat(n,1.0)
    @assert mom.make_mon_expo(n,t,M) == mom.make_mon_expo(n,t)
end









using FileIO

for n in 5:7, k in 2:4
    M = cp_rank.gen_random_band_mat(n,k)
    mc = sort(sparsity.get_maximal_cliques(M))
    mon_cliques = sparsity.get_monomial_cliques(n,t,mc)

    mat_support = sparsity.show_mat_support(repeat(M, inner=(100,100)));
    save("test//support_M_n$(n)_k$(k).png",mat_support)
    support_graph = sparsity.show_support_graph(M);
    save("test//support_graph_n$(n)_k$(k).png",support_graph)
    for t in 2:3
        ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,conlist)
        rmms = cp_model.rec_mom_mat(n,t,ξₜᶜᵖ,mon_cliques)

        mon_mat_support = sparsity.show_mat_support(repeat(rmms, inner=(20,20)));
        save("test//support_mon_mat_n$(n)_k$(k)_t$(t).png",mon_mat_support)
    end
end


n =5
k =2
t =3
conlist = "sGid"
M = cp_rank.gen_random_band_mat(n,k)
mc = sort(sparsity.get_maximal_cliques(M))
mon_cliques = sparsity.get_monomial_cliques(n,t,mc)

mat_support = sparsity.show_mat_support(repeat(M, inner=(100,100)));
save("test//support_M_n$(n)_k$(k).png",mat_support)
support_graph = sparsity.show_support_graph(M);
save("test//support_graph_n$(n)_k$(k).png",support_graph)

ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,conlist)
rmms = cp_model.rec_mom_mat(n,t,ξₜᶜᵖ,mon_cliques)

mon_mat_support = sparsity.show_mat_support(repeat(rmms, inner=(20,20)));
save("test//support_mon_mat_n$(n)_k$(k)_t$(t).png",mon_mat_support)

unique(mom.make_mon_expo(n,(t,t))[(rmms .> 0)])
mom.make_mon_expo(n,2*t)


unique([mom.make_mon_expo(cat(mon_cliques[1:end-1]...,dims=1))...])


MM = mom.make_mon_expo(n,t)
catman  = map(mon -> findall([mon in m for m in mom_cliques]), MM)
sparnar = sortperm(catman)
sMM = mom.make_mon_expo(MM[sparnar])
rec_mom_mat(sMM, ξₜᶜᵖ)



mag = 10
mag_rec_mom_mat = repeat(rec_mom_mat_sparse,inner=(mag,mag))


# show_mat_heat(M)




fat_M = sparsity.show_mat_support(repeat(M,inner=(100,100)));
save("M.png", fat_M)
fat_MOM_3 = sparsity.show_mat_support(mag_rec_mom_mat);
save("MOM_2.png", fat_MOM_2)
fat_MOM_2
save_mat(rec_mom_mat_sparse,"fat_MOM_3.csv")
save("MOM_3.png", fat_MOM_3)

typeof(nar)



ℐₖs_noz = [ℐₖs[i][2:end] for i in 1:3]
ℐₖs_noz[1]
setdiff(ℐₖs_noz[1],ℐₖs_noz[2])

nar = []
for i in 1:2
    push!(nar, setdiff(intersect(ℐₖs_noz[1:i]...) ,union(ℐₖs_noz[(i+1):3]...)))
end
push!(nar, intersect(ℐₖs_noz[1:3]...))
push!(nar, setdiff(intersect(ℐₖs_noz[2:3]...),ℐₖs_noz[1]))
push!(nar, setdiff(ℐₖs_noz[3], union(ℐₖs_noz[1:2]...)))





## running some numerics
n = 6
mats_dir = save_dir*"cp_mats_$n\\"
mkdir(mats_dir)
for k in 1:10
    ext = "Ex_$k"
    rsM_n = cp_matrices.get_random_sparse_cp_mats(n)
    cp_rank.save_mat(rsM_n, mats_dir*ext*"rsM_$n.csv")
end
t = 3
cons = "sGid"
cp_mats = readdir(mats_dir)
comp_dir = save_dir*"\\xi_$cons\\"
mkdir(comp_dir)
for mat in cp_mats
    M = cp_rank.load_mat(mats_dir*mat)
    save_path = comp_dir*mat*"_t$(t)_$cons.txt"
    ξₜᶜᵖ = cp_rank.run_and_save_get_ξₜᶜᵖ(M,t,cons, save_path)
end


ξₜᶜᵖrsM_n_t_sGid = cp_rank.run_and_save_get_ξₜᶜᵖ(rsM_n ,t, save_path)
ξₜᶜᵖ_save_path = save_dir*ext*"ξₜᶜᵖrsM_n$(n)_t$(t)_$cons.txt"

mom_save_path = save_dir*ext*"mom_mat_ξₜᶜᵖrsM_n$(n)_t$(t)_$cons.csv"
save_mom_mat(n,t,ξₜᶜᵖrsM_n_t_sGid, mom_save_path)

# save_mom_mat(n,t,ξₜᶜᵖrsM_n_t_cons)
n = 6; k =3; t = 3
M = cp_matrices.gen_random_band_mat(n,k)




