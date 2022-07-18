using Plots
using LinearAlgebra ; const la = LinearAlgebra

src_dir = dirname(@__FILE__)*"\\"
include(src_dir*"moments.jl")
include(src_dir*"nn_matrices.jl")
include(src_dir*"nn_model.jl")
include(src_dir*"matrix_IO.jl")

include(src_dir*"batch_run.jl")

using .moments ; const mom = moments
using .nn_model ; const nnm = nn_model
#--------------------------------------------------------------------------------------------------
M1 = nn_matrices.make_NN_mat(5)
M2 = nn_matrices.make_lit_NN_mat(1,0.75)
NN_mats = nn_matrices.make_lit_NN_mats()
M3 = NN_mats["Gil1.2"]
M4 = NN_mats["S-MBex4.9"]
M5 = nn_matrices.make_Beasley_Laffey_NN_mat(5)
#--------------------------------------------------------------------------------------------------
M2 = nn_matrices.make_lit_NN_mat(1,0.99)
M3 = nn_matrices.make_lit_NN_mat(1,1)
t = 1


nn_model.computeξₜⁿⁿ(nn_model.modelξₜⁿⁿⁱᵈ(M3, t, true, true, false))
nn_model.computeξₜⁿⁿ(nn_model.modelξₜⁿⁿⁱᵈ(M3, t, true, true, true))

ξₜⁿⁿⁱᵈM  = nn_model.get_ξₜⁿⁿ(M2, t,"id dag ideal")
ξₜⁿⁿⁱᵈM  = nn_model.get_ξₜⁿⁿ(M3, t,"id dag ideal")

ξₜⁿⁿⁱᵈM  = nn_model.get_ξₜⁿⁿ(M2, t,"id ddag")
ξₜⁿⁿⁱᵈM  = nn_model.get_ξₜⁿⁿ(M3, t,"id ddag")
#--------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------

ξₜⁿⁿⁱᵈM  = nn_model.get_ξₜⁿⁿ(M, t,"id")
ξₜⁿⁿⁱᵈdM  = nn_model.get_ξₜⁿⁿ(M, t,"id dag")
ξₜⁿⁿⁱᵈddM = nn_model.get_ξₜⁿⁿ(M,t,"id ddag ")
ξₜⁿⁿˢᵖM = nn_model.get_ξₜⁿⁿ(M,t,"sp")
ξₜⁿⁿˢᵖdM = nn_model.get_ξₜⁿⁿ(M,t,"sp dag")
ξₜⁿⁿˢᵖddM = nn_model.get_ξₜⁿⁿ(M,t,"sp ddag")

ξₜⁿⁿˢᵖM - ξₜⁿⁿⁱᵈM 
ξₜⁿⁿˢᵖdM - ξₜⁿⁿⁱᵈdM 
ξₜⁿⁿˢᵖddM - ξₜⁿⁿⁱᵈddM

#--------------------------------------------------------------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
exa = ["nn"][1]*"\\"
dataexadir = data_dir*exa

# nn_matrices.generate_lit_nn_mats(0.01, dataexadir)
#--------------------------------------------------------------------------------------------------
t = 3
results_dir = assets_dir*"results\\$exa"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_subdir = results_dir*"t$(t)\\"
!isdir(results_subdir) ? mkdir(results_subdir) : 0
mats = matrix_IO.load_mats(dataexadir)
flavs = "nn ddag ".*["id"] # "sp", 


# batch_run.batch_comp_and_save(mats, t, flavs, results_subdir)
batch_run.make_summary(results_subdir)
df = batch_run.clean_summary(results_subdir)

using CSV, DataFrames
df_1 =  CSV.read(results_dir*"t1\\clean_summary.csv", DataFrame,delim =",")
df_2 =  CSV.read(results_dir*"t2\\clean_summary.csv", DataFrame,delim =",")
df_3 =  CSV.read(results_dir*"t3\\clean_summary.csv", DataFrame,delim =",")

bounds_mat =  hcat(Matrix(df_1[:,["nn ddag id", "nn ddag sp"]]),
                Matrix(df_2[:,["nn ddag id", "nn ddag sp"]]),
                Matrix(df_3[:,["nn ddag id", "nn ddag sp"]]))
 
###
disc_blocks = Int.(bounds_mat' .> 1) .+  Int.(bounds_mat' .> 2) .+  Int.(bounds_mat' .> 3) .+  Int.(bounds_mat' .>= 4)            

using  Plots.PlotMeasures 

heatmap(#bounds_mat' ,
        disc_blocks,
        xrotation = 65,
        xticks=(1:5:101,0:0.05:1),
        yticks=(1:6, ["id t=1" "sp t=1" "id t=2" "sp t=2" "id t=3" "sp t=3"],),
        xaxis = "a",
        yaxis = "bound and level",
        title = "ξₜⁿⁿˢᵖ and ξₜⁿⁿⁱᵈ bounds for t = 1,2,3 against a ∈ [0,1] ",
        c =cgrad([:white, :yellow, :red, :purple], [1, 2, 3, 4]),
        # c =cgrad([:white, :black], [0, 4]),
        size = (1000,500),
        bottom_margin = 10mm,
        left_margin = 10mm,
        legend=false)


 #--------------------------------------------------------------------------------------------------

#----- scraps------------
load_dir = results_dir*"t3\\"
list_o_mats = [s for s in readdir(load_dir) if contains(s,".txt")]


check_a(str) = length(split(str,"_")[3][2:end])



for str in list_o_mats
        if check_a(str) < 4
                spl = split(str,"_")
                spl[3] = spl[3]*"0"
                spl = join(spl,"_")
                mv(results_dir*"t3\\"*str, results_dir*"t3\\"*spl)
        end
end

