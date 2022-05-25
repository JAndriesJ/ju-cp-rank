

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"extract_atoms.jl")
include(proj_dir*"cp_matrices.jl")

using .extract_atoms
using .cp_matrices 

datadir = dirname(dirname(proj_dir))*"\\assets\\data\\"
r_cp_mats = cp_matrices.load_mats(datadir*"rand\\")
r_cp_mat = r_cp_mats["ex02_n5_zd0.50_r5_ucpr10_"]

t = 3
flav = "rand"
load_dir = pwd()*"\\assets\\2022-05-10T08\\t$(t)\\$flav\\moments\\"
load_path = load_dir*"Mom_ex02_n5_zd0.50_r5_ucpr10__id_t3.csv"
n = parse(Int,split(file_name,'_')[3][2])
atoms = extract_atoms.get_atoms(n,t,load_path,1e-4)

a = [a.center for a in atoms.atoms ]
w = [a.weight for a in atoms.atoms ]
rec_r_cp_mat = sum([w[k]*a[k]*a[k]'  for k in 1:9])

sum(abs.(r_cp_mat - rec_r_cp_mat))











