using Random
Random.seed!(373)

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"constraints.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"compute.jl")
include(proj_dir*"cp_rank.jl")

using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .constraints ; const con = constraints
using .cp_model
using .compute 
using .cp_rank

using Dates

#generate matrices -------------------------------------------------------------------
Bomze_mats_dict = cp_rank.get_Bomze_cp_mats() # load a dictionary of dense cp matrices from literature.
M = Bomze_mats_dict["M7t"] # pick a matrix from the dictionary
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
M = random_cp_mat
t = 2 # level of the hierarchy
conlist = "sGid" #a string that if it contains sG it impliments the G⊗L ⪰ 0 constriants, If it contains id it impliments the ideal constraints.
ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,conlist)
# Compute ξₜᶜᵖ(M) and save output
save_dir = "C:\\Users\\jandr\\Dropbox\\Vaults\\Pandemonium\\LABOURS\\4-Toulouse\\Numerics\\"*string(now())[1:10]*"\\"
mkdir(save_dir)
cp_rank.run_and_save_get_ξₜᶜᵖ(M,t,conlist,save_dir*"ξₜᶜᵖrsM_n$(n)_t$(t)_$conlist.txt")
## --------------------------------------------



n = 6
t = 3
cons = "sGid"
mats_dir = save_dir*"cp_mats_$n\\"
mkdir(mats_dir)
for k in 1:10
    ext = "Ex_$k"
    rsM_n = cp_matrices.get_random_sparse_cp_mats(n)
    cp_rank.save_mat(rsM_n, mats_dir*ext*"rsM_$n.csv")
end
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

using Graphs, GraphPlot, Fontconfig, Cairo
G = Graphs.Graph(M .> 0)
# gplot(G, nodelabel=[Graphs.vertices(G)...])
MC = Graphs.maximal_cliques(G)
MC = sort(map(x ->sort(x),MC))

ℐₖs = get_ℐₖs(t,MC,n)

get_ℐₖs(t,cc,n) = [ get_ℐₖ(t,cc[k],n) for k in 1:length(cc)]
get_ℐₖ(t,c,n) = map(v->embed(v,c,n), make_mon_expo(length(c),t))
embed(v,α,n) = [i ∈ α ? popfirst!(v) : 0  for i in 1:n]
get_ℐₖ_compliment(mon, ℐₖ) =  setdiff(mon, ℐₖ)


ℐₖ = unique(cat(ℐₖs...,dims=1))
mon_ex = make_mon_expo(n,t)
ℐₖ_compliment = get_ℐₖ_compliment(mon_ex, ℐₖ)

nar = cat(ℐₖ,ℐₖ_compliment,dims=1)


A = [α+β for α ∈ nar, β ∈ nar]

ξₜᶜᵖ = cp_rank.get_ξₜᶜᵖ(M,t,"sGid")
MOM = cp_model.rec_mom_mat(A,ξₜᶜᵖ)

B = Int.(MOM .> 0)
B = repeat(B,inner = (10,10))

pan = Gray.(B)
pan = Gray.(Int.(M .> 0))

# @assert Int.(sum([ Xᵢⱼ(n,e[1],e[2]) for e in nze]) + I(6) .> 0) == M
#-------------------------------------------------------------------
# include(proj_dir*"Batch_proc.jl")
# using .Utility


cp_mats = cpMatrices.get_cp_mats()  # load the matrices
A = cp_mats["M6"]                   # Pick a specific one
t = 2                               # Choose the level of the hierarchy
conlist = "DagXXwGsG"               # Choose which additional constraints to add
mod = cpModel.Modelξₜᶜᵖ(A,t,conlist) # Build the model
Compute.Computeξₜᶜᵖ(mod)  # Solve the model using Mosek

