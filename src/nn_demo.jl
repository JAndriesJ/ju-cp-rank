#-------------------------------------Diagnostics--------------------------------------------
test_dir = dirname(dirname(@__FILE__))*"\\test\\"
include(test_dir*"nn_matrices_test.jl")
include(test_dir*"nn_model_test.jl")

#-------------------------------------Initialization--------------------------------------------
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
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
using JuMP
t = 1
# M_nar(n) = [i ≥ j ? 1 : 0  for i ∈ 1:n, j ∈ 1:n ] .+ 0.0
M_nar(n) =  diagm(ones(n)) .+ 0.0
function get_m_ext(n)
    M =  M_nar(n)#nn_matrices.make_nn_ED_mat(n) ; 
    _, ext_mom    = nn_model.get_ξₜⁿⁿ(M, t,"id") ; 
    ext_mom_dic = Dict(zip(ext_mom[:,1],ext_mom[:,2]))
    M_ext = map(m -> ext_mom_dic[m],  moments.make_mon_expo(moments.make_mon_expo(1, M, true)))
    return M_ext #round.(M_ext, digits=3)
end
get_m_ext(4)

m_ext_mod = [ get_m_ext(n+4) for n in 1:3]
n = 10
A = deepcopy( get_m_ext(n) )
M = A[2:n+1,n+2:2n+1]

C = cholesky(A[2:end,2:end])
U1 = C.U[1:10,1:10] #.* (C.U .> 1.0e-6)
U1'*U1
A

U^2
N = U'*U
A - N


N = hcat(vcat(2*diagm([16,9,4,9,16]),M),vcat(M,2*diagm([16,9,4,9,16]))  )
x = diag(N)/(n-1)
B = vcat([8 x'], hcat(x,N))


M = [1 1 0 0 0 0 1
     1 1 1 0 0 1 0
     0 1 1 1 1 0 0
     0 0 1 1 1 0 0
     0 0 0 1 1 1 0
     0 0 0 0 1 1 1
     0 0 0 0 0 1 1]

M = [0 1 0 0 0 0 0
     1 0 1 0 0 0 0
     0 1 0 1 0 0 0
     1 0 1 0 1 0 0
     0 1 0 1 0 1 0
     1 0 1 0 1 0 1
     0 1 0 1 0 1 0] .+ 0.


M = [1 0 0 0 0 0 0 0
     0 1 0 0 0 0 0 0
     1 0 1 0 0 0 0 0
     0 1 0 1 0 0 0 0
     1 0 1 0 1 0 0 0
     0 1 0 1 0 1 0 0
     1 0 1 0 1 0 1 0
     0 1 0 1 0 1 0 1] .+ 0.



M = [1 1 0 1 0 1 0 1
     0 1 1 0 1 0 1 0
     1 0 1 1 0 1 0 1
     0 1 0 1 1 0 1 0
     1 0 1 0 1 1 0 1
     0 1 0 1 0 1 1 0
     1 0 1 0 1 0 1 1
     0 1 0 1 0 1 0 1] .+ 0.

M_nar(n) = [i ≥ j ? 1 : 0  for i ∈ 1:n, j ∈ 1:n ] .+ 0.0
n = 10
M = M_nar(n)
J = ones(n,n)
I = 3 .* diagm(ones(n))

rank(M)


A = vcat(hcat(J+I,M), hcat(M',J+I))
eigvals(A)

B = vcat(hcat(n,3*ones(2n)'),hcat(3*ones(2n),A))
eigvals(B)

# tests
# N[2:n+1,n+2:2n+1] #(57)
# eigvals(A) # (58)
# all([(n-1)*N[1,i] - N[i,i] for i in 2:2n+1] .≥ 0) # (59)
# [ MK_ultra[i,j]*MK_ultra[1,1] - MK_ultra[i,j] for i in 2:n+1, j ∈ n+2:2n+1] # (60)




# ξₜⁿⁿⁱᵈdM, _   = nn_model.get_ξₜⁿⁿ(M, t,"id dag")
ξₜⁿⁿⁱᵈddM, _  = nn_model.get_ξₜⁿⁿ(M,t,"id ddag")
# ξₜⁿⁿˢᵖM, _    = nn_model.get_ξₜⁿⁿ(M,t,"sp")
# ξₜⁿⁿˢᵖdM, _   = nn_model.get_ξₜⁿⁿ(M,t,"sp dag")
ξₜⁿⁿˢᵖddM, _  = nn_model.get_ξₜⁿⁿ(M,t,"sp ddag")

ξₜⁿⁿˢᵖM - ξₜⁿⁿⁱᵈM 
ξₜⁿⁿˢᵖdM - ξₜⁿⁿⁱᵈdM 
ξₜⁿⁿˢᵖddM - ξₜⁿⁿⁱᵈddM


#--------------------------------------------------------------------------------------------------
assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\"
exa = ["nn" "EDM"][2]*"\\"
dataexadir = data_dir*exa
# nn_matrices.generate_lit_nn_mats(0.01, dataexadir)
# nn_matrices.generate_ED_mats([4:20...], dataexadir)


M = nn_matrices.make_nn_ED_mat(32)
ξₜⁿⁿⁱᵈdM, _   = nn_model.get_ξₜⁿⁿ(M, 1,"id dag")
#--------------------------------------------------------------------------------------------------
mats = matrix_IO.load_mats(dataexadir)
t = 3
flavs = "nn ddag ".*["id", "sp"] # , 

results_dir = assets_dir*"results\\$exa"
!isdir(results_dir) ? mkdir(results_dir) : 0
results_subdir = results_dir*"t$(t)\\"
!isdir(results_subdir) ? mkdir(results_subdir) : 0

batch_run.batch_comp_and_save(mats, t, flavs, results_subdir)
batch_run.make_summary(results_subdir)
df = batch_run.clean_summary(results_subdir)
#--------------------------------------------------------------------------------------------------
batch_run.get_nn_ex_fig(results_dir)



