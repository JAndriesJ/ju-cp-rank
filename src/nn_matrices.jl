module nn_matrices
using Random

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")

export make_NN_mat,
       make_lit_NN_mat,
       make_lit_NN_mats,
       make_Beasley_Laffey_NN_mat
    #    make_NN_supp_mat


# function make_NN_supp_mat(r,c,q = rand(1:n))
#     n = r*c
#     zov = shuffle(vcat(ones(Int,q), zeros(Int,n-q)))
#     return reshape(zov,(r,c))
# end
make_Beasley_Laffey_NN_mat(n) = [(i-j)^2 for i ∈ 1:n, j∈ 1:n] .+ 0.0


make_NN_mat(n) = [ (i==j+1 || i==j-1 || i==j==1 || i==j==n ) ? 1.0 : 0.0 for i ∈ 1:n, j ∈ 1:n ] .+ 0.0
function make_NN_mat(r,c,q = rand(1:r*c))
    zov = round.(shuffle(vcat(rand(q).+0.1, zeros(r*c-q))),digits=2)
    return reshape(zov,(r,c)) 
end
make_NN_mat() = make_NN_mat(rand(3:8,2)...)


function make_lit_NN_mats()
    nn_dict = Dict()
    nn_dict["Gil1.2"] = [1 1 0 0
                         1 0 1 0
                         0 1 0 1
                         0 0 1 1 ] .+ 0.0

    nn_dict["S-MBex4.9"] = [6 3 3 0
                            3 5 1 3
                            3 1 5 3
                            0 3 3 6] .+ 0.0

    # nn_dict["Gil2.6"] = make_NN_hex() # dense
    return nn_dict
end

"""
Lower Bounds on Matrix Factorization Ranks via
Noncommutative Polynomial Optimization
Sander Gribling1 · David de Laat1 · Monique Laurent1,2
"""
make_lit_NN_mat(a,b) =     [1-a 1+a 1-b 1+b
                            1+a 1-a 1-b 1+b
                            1+a 1-a 1+b 1-b
                            1-a 1+a 1+b 1-b] .+ 0.0
function  generate_lit_nn_mats(step_size::Float64,save_dir::String)
    for a ∈ 0:step_size:1.0
        b = 1
        mat = make_lit_NN_mat(a,b)
        a = rpad(a, 4, '0')
        matrix_IO.save_mat(mat, save_dir*"ex_a$(a)_b$(b)_.csv")
    end
end


function  generate_ED_mats(n_range,save_dir::String)
    for n ∈ n_range
        mat = make_nn_ED_mat(n)
        matrix_IO.save_mat(mat, save_dir*"ex_EDM_n$(n)_.csv")
    end
end
make_nn_ED_mat(n) = [ (i-j)^2 for i ∈ 1:n, j ∈ 1:n] .+ 0.0

make_NN_hex(a=2) =   circmatrix([1, a, 2a-1,  2a-1,   a,  1])
circmatrix(v) = hcat([circshift(v, k) for k = -length(v)+1:0]...)



end


