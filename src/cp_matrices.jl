module cp_matrices
using Test
using LinearAlgebra
using Random

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")
using .matrix_IO 

export  get_lit_cp_mats,
        get_random_cp_mats,
        get_random_sparse_cp_mats,
        gen_random_band_mat,
        run_tests

## Examples

function get_lit_cp_mats()
    lit_cp_mats = Dict()
    #  Naomi Shaked-Monderer and Abraham Berman. “Copositive and Completely Positive Matrices.” (2021): 265-302
    lit_cp_mats["ex01"] = [1 1 0 0 1
                      1 2 1 0 0
                      0 1 3 1 0
                      0 0 1 4 1
                      1 0 0 1 5]

    lit_cp_mats["ex02"] = [2 1 1 0 0
                      1 2 0 0 1
                      1 0 2 1 0
                      0 0 1 2 1
                      0 1 0 1 2]

    # Francesco Barioli. “Completely positive matrices with a book-graph.”  Linear Algebra and its Applications 277 (1998): 11-31  
    lit_cp_mats["ex03"] = [3 2  0 0 1
                      2 5  6 0 0
                      0 6 14 4 0
                      0 0  4 9 1
                      1 0  0 1 2]

    lit_cp_mats["ex04"] =[12  2 2  6 4 2
                      2 23 7  9 6 4
                      2  7 5  7 0 0
                      6  9 7 13 0 0
                      4  6 0  0 8 6
                      2  4 0  0 6 5]

    lit_cp_mats["ex05"] =[13  0 1 2 1 1 4 1 1 1
                      0 19 3 3 1 3 1 1 3 5
                      1  3 5 4 0 0 0 0 0 0
                      2  3 4 5 0 0 0 0 0 0
                      1  1 0 0 2 3 0 0 0 0
                      1  3 0 0 3 5 0 0 0 0
                      4  1 0 0 0 0 5 2 0 0
                      1  1 0 0 0 0 2 1 0 0
                      1  3 0 0 0 0 0 0 2 4
                      1  5 0 0 0 0 0 0 4 10]
      
    # Turkel A. Haufmann. “On Completely Positive Matrices.” Master’s Thesis for the degree of Master of Applied Mathematics and Mechanics (2011): 30-31
    lit_cp_mats["ex06"] = [2 0 1  5  0  0
                      0 1 1  0  0  0
                      1 1 2  0  0  0
                      5 0 0 31  5  6
                      0 0 0  5 13  8
                      0 0 0  6  8 17]

    # Luigi Salce and Paolo Zanardo “Completely Positive Matrices and Positivity of Least Squares Solutions.” (1993)  
    lit_cp_mats["ex07"] =  [ 2 1 0 0 1
                        1 2 1 0 0
                        0 1 2 1 0
                        0 0 1 2 1
                        1 0 0 1 3]

    #  Shuhuang Xiang Xi'an, and Shuwen Xiang. “Notes on Completely Positive Matrices.” Linear Algebra and its Applications 271 (1998) 273-282  
    lit_cp_mats["ex08"] = [2 0 0 1 1
                      0 2 0 1 1
                      0 0 2 1 1
                      1 1 1 3 0
                      1 1 1 0 3]

    lit_cp_mats["ex09"] =[15   3  2 1 0
                      3  18  7 0 2
                      2   7 20 6 1
                      1   0  6 9 2
                      0   2  1 2 7]

    lit_cp_mats["ex10"] =[2 0 1 1 0
                      0 2 1 1 0
                      1 1 4 1 0
                      1 1 1 5 0
                      0 0 0 0 3]  

    #  Abraham Berman & Naomi Shaked-Monderer. “Completely Positive Matrices: Real, Rational, and Integral.” (2017) 
    lit_cp_mats["ex11"] =[15  3  2  1  0
                       3 18  7  0  2
                       2  7 20  6  1
                       1  0  6  9  2
                       0  2  1  2  7]
      
    # Hohgbo Dong & Kurt Anstreicher. “A note on 5 x 5 Completely Positive Matrices” Linear Algebra and its Applications 433 (2010) 1001-1004 
    lit_cp_mats["ex12"] = [8 5 1 1 5
                       5 8 5 1 1
                       1 5 8 5 1
                       1 1 5 8 5
                       5 1 1 5 8]

    lit_cp_mats["ex13"] = [2 0 0 1 0
                       0 2 0 0 1
                       0 0 4 1 1
                       1 0 1 1 0
                       0 1 1 0 1]
    
    # Abraham Berman & Changqing Xu. “5 x 5 Completely positive matrices with a book-graph.”  Linear Algebra and its Applications 393 (2004): 55-71
    lit_cp_mats["ex14"] = [6 5 5 1 2
                       5 6 5 2 1
                       5 5 9 2 2
                       1 2 2 1 0
                       2 1 2 0 1]

    # Mikolo Cedolin & Luigi Salce. “Completely positive matrices of order 5 with CP-graph.” Note Mat. 36 (2016) no. 1, 123-132
    lit_cp_mats["ex15"] =  [10  5  2  0  5
                         5 10  5  0  0
                         2  5 10  5  0
                         0  0  5 10  5
                         5  0  0  5 10]
   
    lit_cp_mats["ex16"] =[10  5   2  0    5
                       5 10   5 10/3  0
                       2  5  10  5    0
                       0 10/3 5 10    5
                       5  0   0  5   10]
        
    # Radu loan Bot & Dang-Khoa Nguyen. “Factorization of completely positive matrices using iterative projected gradient steps.” (2021) and P. Groetzner, M. Dur, A factorization method for completely positive matrices, Linear Algebra Appl. (2020).
    lit_cp_mats["ex17"] =[41 43  80  56 50
                      43 62  89  78 51
                      80 89 162 120 93
                      56 78 120 104 62
                      50 51  93  62 65] 

     # Burer, Samuel; Antreicher, Kurt M.; Duer, Mirjam. “The difference between 5 x 5 doubly nonnegative and completely positive matrices.” Linear Algebra and its Applications 431(9) (2009) 1539-1552
    lit_cp_mats["ex18"] = [24 11   0   0 11
                       11 24  11   0  0
                        0 11 124  11  0
                        0  0  11  24 11
                       11  0   0  11 24]

    # Max Pfeffer & Jose Alejandro Samper. “The cone of 5 x 5 Completely Positive Matrices.”  (2021)
    lit_cp_mats["ex19"] = [53 32  1  4  26
                       32 30  5  2   5
                        1  5 11 10   3
                        4  2 10 21  17
                       26  5  3 17  35]

    # 2014, Nie, The A-Truncated K-Moment Problem
    lit_cp_mats["ex20"] =[6 4 1 2 2
                      4 6 0 1 3
                      1 0 3 1 2
                      2 1 1 2 1
                      2 3 2 1 5] 
    # from Completely positive matrices, Berman Shaked-Monderer,  World scientific
    lit_cp_mats["ex21"] = [  8  12 16 4 6  8
                      12 20 28 6 10 14
                      16 28 40 8 14 20
                      4   6  8 2  3  4
                      6  10 14 3  5  7
                      8  14 20 4  7 10]
    #
    lit_cp_mats["ex22"] = [13.7162 8.1090 10.0752  7.3940 10.7332 4.2551 9.8380
                        8.1090 9.9804  6.7089  6.1420  8.3044 3.9288 6.7772
                       10.0752 6.7089 69.6070  5.6133  8.3782 3.6309 6.6191
                        7.3940 6.1420  5.6133 10.2718  8.1573 4.3024 5.9480
                       10.7332 8.3044  8.3782  8.1573 12.8697 4.1992 9.1688
                        4.2551 3.9286  3.6309  4.3024  4.1092 2.5809 3.3640
                        9.8380 6.7772  6.6191  5.9480  9.1688 3.3640 13.1836]
    
    
    # Matrices from:  I.M. Bomze, W. Schachinger, R. Ullrich. From seven to eleven: Completely positive matrices with high cp-rank. Linear Algebra and its Applications 459 (2014), 208 – 221.
    lit_cp_mats["ex23"] = circmatrix(([531.0, 81, 150, 150, 81, 531, 926]))
    lit_cp_mats["ex24"] =  circmatrix(([108, 27, 4, 4, 27, 108, 163.0]))
    lit_cp_mats["ex25"] =  [541.0 880 363 24 55 11 24 0
                          880 2007 1496 363 48 22 22 24
                          363 1496 2223 1452 363 24 22 11
                          24 363 1452 2325 1584 363 48 55
                          55 48 363 1584 2325 1452 363 24
                          11 22 24 363 1452 2223 1496 363
                          24 22 22 48 363 1496 2007 880
                          0 24 11 55 24 363 880 541]

    lit_cp_mats["ex26"] =  [2548 1628 363 60 55 55 60 363 1628
                      1628 2548 1628 363 60 55 55 60 363
                      363 1628 2483 1562 363 42 22 55 60
                      60 363 1562 2476 1628 363 42 55 55
                      55 60 363 1628 2548 1628 363 60 55
                      55 55 42 363 1628 2476 1562 363 60
                      60 55 22 42 363 1562 2483 1628 363
                      363 60 55 55 60 363 1628 2548 1628
                      1628 363 60 55 55 60 363 1628 2548]

    lit_cp_mats["ex27"] =  [781 0 72 36 228 320 240 228 36 96 0;
                        0 845 0 96 36 228 320 320 228 36 96
                      72 0 827 0 72 36 198 320 320 198 36
                      36 96 0 845 0 96 36 228 320 320 228
                      228 36 72 0 781 0 96 36 228 240 320
                      320 228 36 96 0 845 0 96 36 228 320
                      240 320 198 36 96 0 745 0 96 36 228
                      228 320 320 228 36 96 0 845 0 96 36
                      36 228 320 320 228 36 96 0 845 0 96
                      96 36 198 320 240 228 36 96 0 745 0
                      0 96 36 228 320 320 228 36 96 0 845]

    for k in keys(lit_cp_mats)
      lit_cp_mats[k] = makediagone(lit_cp_mats[k])
    end
    return lit_cp_mats
end

function get_lit_non_cp_mats()
  lit_non_cp_mats = Dict()
## 2014, Nie, The A-Truncated K-Moment Problem
  lit_non_cp_mats[1] = [1 1 0 0 1
                        1 2 1 0 0
                        0 1 2 1 0
                        0 0 1 2 1
                        1 0 0 1 6] 
  for k in keys(lit_cp_mats)
    lit_cp_mats[k] = makediagone(lit_cp_mats[k])
  end                      
  return lit_cp_mats
end

function get_random_cp_mats(n::Int,r::Int)
    a = rand(n,r)
    return r, a, makediagone(sum([a[:,k]*a[:,k]' for k ∈ 1:r ]))   
end

function get_random_cp_mats(M::Matrix{Int64})
    n = size(M)[1]
    nze = get_nonzero_entries(M)
    R = cat([ reᵢⱼ(n,e[1],e[2]) for e in nze]...,dims=2)
    r = size(R)[2]
    return r,R, makediagone(R*R' .+ 0.0)
end
get_random_cp_mats(n::Int,p::Float64=0.5) = get_random_cp_mats(gen_random_sparcity_mat(n,p))
gen_random_sparcity_mat(n,p=0.5) = Int.(Symmetric(map(x -> bernoulli(p), ones(n,n))) + I(n) .> 0)
reᵢⱼ(n::Int,i::Int,j::Int) = [k ∈ [i,j] ? rand() : 0 for k in 1:n]
# rXᵢⱼ(n::Int,i::Int,j::Int) = reᵢⱼ(n,i,j)*reᵢⱼ(n,i,j)'
get_nonzero_entries(M) = [(i,j) for i in 1:(size(M)[1]-1) for j in i:size(M)[1] if M[i,j] != 0]
bernoulli(p) =  rand() <= p ? 1 : 0

# Generates a batch of cp-matrices, stores them, their construciton, and rank 
function get_random_cp_mats(n_range, num_ex, save_dir)
  Random.seed!(373)
    for n in n_range, 
        p = rand()/2 +0.25
        for ex in num_ex
            ex = string(ex,pad=2)
            rp,R,M = get_random_cp_mats(n,p)
            r = LinearAlgebra.rank(M .+ 0.0)
            M_name = "M_n$(n)_ex$(ex)_r$(r)_rp$(rp).csv"
            matrix_IO.save_mat(M, save_dir*M_name)
            R_name = "R_n$(n)_ex$(ex)_r$(r)_rp$(rp).csv"
            matrix_IO.save_mat(R, save_dir*R_name)
        end
    end
end
# Loads all the matrices stored as .csv's from a directory as a dictionary
function get_random_cp_mats(load_dir::String)
    list_o_mats = readdir(load_dir)
    list_o_mat_names = [m[1:end-4] for m in list_o_mats]
    mats = [matrix_IO.load_mat(load_dir*mat) for mat in list_o_mats]
    return Dict(zip(list_o_mat_names, mats))
end

function get_random_NN_band_mat(n::Int, k::Int, r=rand(3:8))
    A = zeros(n,n)
    for s ∈ 0:(n-k) 
        B = cat(zeros(s,r),rand(k,r),zeros(n-k-s,r),dims=1)
        A += B*B'
    end
    return A
end
function get_random_NN_band_mat(n_range, k_range)
  @assert length(n_range) == length(k_range)
  @assert n_range .> k_range
  return [get_random_NN_band_mat(n, k) for n in n_range, k in k_range ]
end

### Utility
"""Vector input, output a Circulant matrix of said vector"""
circmatrix(v) = hcat([circshift(v, k) for k = -length(v)+1:0]...)
"""Makes the diagonal all ones via scaling"""
makediagone(A) = diagm(1 ./ sqrt.(diag(A))) * A * diagm( 1 ./ sqrt.(diag(A)))

isPSD(A) = minimum(LinearAlgebra.eigvals(A)) > -2.e-14
isNN(A)  = all(A .≥ 0)
isDNN(A) = isPSD(A) && isNN(A)
isDiagOne(A) = -1e-15 ≤  sum(diag(A))-size(A)[1] ≤ 1e-15

function CP_necessities(A) 
  @test isPSD(A)
  @test isNN(A) 
  @test isDiagOne(A)
end 
  
function run_tests()
  @testset "Misc." begin
    @test circmatrix([1,2,3])  == [3 2 1; 1 3 2; 2 1 3]
    @test makediagone([1,2,3,4]*[1 2 3 4]) == [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
  end

  @testset "get_lit_cp_mats" begin 
    lit_cp_mats = get_lit_cp_mats()  
    @test length(lit_cp_mats) == 27
    for k in keys(lit_cp_mats)
        CP_necessities(lit_cp_mats[k]) 
    end
  end

  @testset "load_random_cp_mats" begin
    load_dir = dirname(dirname(@__FILE__))*"\\assets\\data\\randomly_generated_cp-matrices\\mats\\"
    random_cp_mats = load_random_cp_mats(load_dir)
    @test length(random_cp_mats) == 50
    for k in keys(random_cp_mats)
      CP_necessities(random_cp_mats[k]) 
    end
  end

  # @testset "get_random_cp_mats(n,r)" begin
  #     for n in 5:10, r in rand(1:20,5)
  #         r,R,A = get_random_cp_mats(n,r)
  #         @test isNN(A) 
  #         @test size(A)[1] == n
  #         @test isPSD(A)
  #     end
  # end    

end

end 
