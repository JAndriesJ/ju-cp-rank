module Utility
using LinearAlgebra
using JuMP

    export  circmatrix,
            makediagone,
            var_kron,
            eₖ,
            index_to_var

## Examples
"""Vector input, output a Circulant matrix of said vector"""
function circmatrix(v)
    hcat([circshift(v, k) for k = -length(v)+1:0]...)
end

"""Makes the diagonal all ones via scaling"""
makediagone(A) = diagm(1 ./ sqrt.(diag(A))) * A * diagm( 1 ./ sqrt.(diag(A)))

## Moments
"""
input: Dictionary:
    keys: i,j ∈ [n], n::Int
    values: square array A_{i,j}
output: Array A = (A_{i,j})_i,j
"""
function assemble_dict(dict_of_blocks)
    n = Int(sqrt(length(keys(dict_of_blocks))))
    if n == 1
        return dict_of_blocks[1,1]
    end
    block = []
    row_block = []
    for i in 1:n
        for j in 1:n
            if j == 1
                row_block = dict_of_blocks[i, j]
            else
                row_block = hcat(row_block, dict_of_blocks[i,j])
            end
        end

        if i == 1
            block = row_block
        elseif i > 1
            block = vcat(block, row_block)
        end

    end
    return block
end

"""
input: A,B (arrays of integer tupples)
output: exponent array of A ⊗ B
"""
function var_kron(A,B)
    n1,n2 = size(A)

    D = Dict()
    for i in 1:n1
        for j in 1:n2
            C = repeat( [A[i,j]] , inner = (1,1), outer = size(B))
            D[i,j] = C + B
        end
    end
    return assemble_dict(D)
end


"""
input: n (integer),k(integer)
output: eₖ ∈ {0,1}ⁿ i.e. the standard basis vector
"""
eₖ(n::Int,k::Int) = vcat(zeros(Int32,k-1),[1],zeros(Int32,n-k))


## Constraints
"""Takes an array of exponents α's and gives array of same shape L(xᵅ)
Input: Array B of multi-indices: α
Output: L of x of array of indices: L(x.(B))
 """
function index_to_var(var, index_array)
    sub = α -> var[α]
    var_array = sub.(index_array)
    return var_array
end

## Model
"""
import the JuMP model to SDPA format
"""

function read_dat_s_model(file_name::String)
    return JuMP.read_from_file(file_name, format=MOI.FileFormats.FORMAT_SDPA)
end

"""
Exmports the JuMP model to SDPA format
"""
function export_model(model,file_name::String)
    JuMP.write_to_file(model, file_name, format=MOI.FileFormats.FORMAT_SDPA)
end



end
