
############################################### clean this shit up Goddammit ########################################
proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"extract_atoms.jl")

assets_dir = dirname(dirname(proj_dir))*"\\assets\\"
mom_dir = assets_dir*"results\\lit\\t3\\moments\\"
mom_dir_list = [c for c in readdir(mom_dir)]

mom_path = mom_dir*mom_dir_list[8]


cen_id, wei_id, A =  extract_atoms.get_atoms(mom_dir*mom_dir_list[7])


x0 = hcat(cen_id.*wei_id...)
x0[x0  .< 1e-6] .= 0
x0 = reshape(x0,1,:)


cen, wei, A =  extract_atoms.get_atoms(mom_path)
A_rec = extract_atoms.recon_mat(cen, wei)


A
A_rec
maximum(abs.(A - A_rec))
sum((A - A_rec).^2)


function read_atoms(atom_path::String)
        df_r = CSV.read(atom_path, DataFrame)
        cen = eval.(Meta.parse.(df_r.centers))
        wei = eval.(Meta.parse.(df_r.weights))
        return cen, wei
end

cen, wei = read_atoms("atoms.csv")
###
n = size(A)[1]
[ (i,j) for i ∈ 1:n for j ∈ i:n if A[i,j] != 0.0 ]


function f_A(A,x0)
        n = size(A)[1]
        p = div(length(x0),n)
        funk(x) = sum([(A[i,j] - sum( [ x[i+k*n]*x[j+k*n] for k ∈ 0:(p-1)  ]))^2  for i ∈ 1:n for j ∈ i:n])
        return funk
end

f = f_A(A,x0)
reshape(x0,5,:)

f(x0)



##### Refinemenet
A =    [2.0  1.0  1.0  1.0  1.0
        1.0  1.0  1.0  1.0  1.0
        1.0  1.0  1.0  1.0  1.0
        1.0  1.0  1.0  1.0  1.0
        1.0  1.0  1.0  1.0  1.0]
#
x0 = vcat(ones(Float64,5), [1, 0.03, 0.03, 0.03, 0.03])- 0.02*rand(10)
p = length(x0)/
function f_A(A,p)
        n = size(A)[1]
        funk(x) = sum([(A[i,j] - sum( [ x[i+k*n]*x[j+k*n] for k ∈ 0:(p-1) ]))^2  for i ∈ 1:n for j ∈ i:n ])
        return funk
end

#--------

f_eval(y) = substitute.(f, (Dict([ (x[i,k],y[i,k] ) for i ∈ 1:n for k ∈ 1:p]),))[1]
function g_eval(y) 
        for z ∈ 1:length(g)
                g_e[z] = substitute.(g[z], (Dict([ (x[i,k],y[i,k] ) for i ∈ 1:n for k ∈ 1:p]),))[1] 
        end
        return g_e
end
H_eval(y) = [substitute.(H[z,a], (Dict([ (x[i,k],y[i,k] ) for i ∈ 1:n for k ∈ 1:p]),))[1] for z ∈ 1:length(g), a ∈ 1:length(g)]

f_eval(x0)
g_eval(x0)
H_eval(x0)

using Optim

Optim.optimize(f_eval, x0)





