module cp_matrices
using Test
using LinearAlgebra
using CSV, DataFrames

export  get_Bomze_cp_mats,
        get_random_cp_mats,
        get_random_sparse_cp_mats,
        gen_random_band_mat,
        get_zero_entries,
        run_tests

## Examples
"""Vector input, output a Circulant matrix of said vector"""
circmatrix(v) = hcat([circshift(v, k) for k = -length(v)+1:0]...)

"""Makes the diagonal all ones via scaling"""
makediagone(A) = diagm(1 ./ sqrt.(diag(A))) * A * diagm( 1 ./ sqrt.(diag(A)))

function get_Bomze_cp_mats()
    cpmatDict = Dict()
    # from Completely positive matrices, Berman Shaked-Monderer,  World scientific
    cpmatDict["M6"] = makediagone(          ([   8  12 16 4 6   8;
                                                12 20 28 6 10 14;
                                                16 28 40 8 14 20;
                                                4   6  8 2  3  4;
                                                6  10 14 3  5  7;
                                                8  14 20 4  7 10]))


    # Matrices from:  I.M. Bomze, W. Schachinger, R. Ullrich. From seven to eleven: Completely positive matrices with high cp-rank. Linear Algebra and its Applications 459 (2014), 208 – 221.
    cpmatDict["M7"] = circmatrix(           ([531.0, 81, 150, 150, 81, 531, 926])/926.0)

    cpmatDict["M7t"] =  circmatrix(         ([108, 27, 4, 4, 27, 108, 163.0])/163.0)

    cpmatDict["M8t"] =  makediagone(         ([541.0 880 363 24 55 11 24 0;
                                               880 2007 1496 363 48 22 22 24;
                                               363 1496 2223 1452 363 24 22 11;
                                               24 363 1452 2325 1584 363 48 55;
                                               55 48 363 1584 2325 1452 363 24;
                                               11 22 24 363 1452 2223 1496 363;
                                               24 22 22 48 363 1496 2007 880;
                                               0 24 11 55 24 363 880 541]))

    cpmatDict["M9t"] =  makediagone(         ([2548 1628 363 60 55 55 60 363 1628;
                                               1628 2548 1628 363 60 55 55 60 363;
                                               363 1628 2483 1562 363 42 22 55 60;
                                               60 363 1562 2476 1628 363 42 55 55;
                                               55 60 363 1628 2548 1628 363 60 55;
                                               55 55 42 363 1628 2476 1562 363 60;
                                               60 55 22 42 363 1562 2483 1628 363;
                                               363 60 55 55 60 363 1628 2548 1628;
                                               1628 363 60 55 55 60 363 1628 2548]))

    cpmatDict["M11t"] =  makediagone(        ([781 0 72 36 228 320 240 228 36 96 0;
                                                 0 845 0 96 36 228 320 320 228 36 96;
                                                72 0 827 0 72 36 198 320 320 198 36;
                                                36 96 0 845 0 96 36 228 320 320 228;
                                                228 36 72 0 781 0 96 36 228 240 320;
                                                320 228 36 96 0 845 0 96 36 228 320;
                                                240 320 198 36 96 0 745 0 96 36 228;
                                                228 320 320 228 36 96 0 845 0 96 36;
                                                36 228 320 320 228 36 96 0 845 0 96;
                                                96 36 198 320 240 228 36 96 0 745 0;
                                                0 96 36 228 320 320 228 36 96 0 845]))
    return cpmatDict
end

function get_random_cp_mats(n,r)
    a = rand(n,r)
    return sum([a[:,k]*a[:,k]' for k ∈ 1:r ])    
end

function get_random_sparse_cp_mats(M::Matrix{Int64})
    n = size(M)[1]
    nze = get_nonzero_entries(M)
    sum([ rXᵢⱼ(n,e[1],e[2]) for e in nze]) + diagm(rand(n))
end
get_random_sparse_cp_mats(n::Int,p::Float64=0.5) = get_random_sparse_cp_mats(gen_random_sparcity_mat(n,p))

reᵢⱼ(n::Int,i::Int,j::Int) = [k ∈ [i,j] ? rand() : 0 for k in 1:n]
rXᵢⱼ(n::Int,i::Int,j::Int) = reᵢⱼ(n,i,j)*reᵢⱼ(n,i,j)'
get_nonzero_entries(M) = [(i,j) for i in 1:(size(M)[1]-1) for j in i:size(M)[1] if M[i,j] != 0]

gen_random_sparcity_mat(n,p=0.5) = Int.(Symmetric(map(x -> bernoulli(p), ones(n,n))) + I(n) .> 0)
bernoulli(p) =  rand() <= p ? 1 : 0


function gen_random_band_mat(n,k)
    @assert k < 5
    B = rand(n,rand(5:10))
    A = B*B'
    A = [abs(i-j) < k ? A[i,j] : 0 for i in 1:n, j in 1:n]
    return A + abs(minimum(eigvals(A)))*I(n)
end

###
isPSD(A) = minimum(LinearAlgebra.eigvals(A)) > -2.e-14

function run_tests()
    @testset "get_Bomze_cp_mats" begin
        Bomze_cp_mats = get_Bomze_cp_mats()
        for mat in Bomze_cp_mats
            @test all(mat .> 0)
            # @test size(A) = (n,)
            @test isPSD(A)
        end
    end


    @testset "get_random_cp_mats(n,r)" begin
        for n in 5:10, r in rand(1:20,5)
            A = get_random_cp_mats(n,r)
            @test all(A .> 0)
            @test size(A) = (n,)
            @test isPSD(A)
        end
    end    


end

end  # module cpMatrices
