module cp_matrices
using Test
using LinearAlgebra
using Random

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")
include(dirname(@__FILE__)*"\\"*"moments.jl")

using .matrix_IO 
using .moments 

export  generate_lit_cp_mats,
        generate_random_cp_mats,
        gen_bipart_supp_mats,
        gen_asc_nzd_cp_mats,
        generate_lit_non_cp_mats,
        load_mats
        # gen_random_band_mat


# Loads all the matrices stored as .csv's from a directory as a dictionary
function load_mats(load_dir::String)
  if load_dir[end-3:end] == ".csv"
    return matrix_IO.load_mat(load_dir) 
  else
    list_o_mats = [s for s in readdir(load_dir) if (contains(s,".csv") && contains(s,"ex"))]
    list_o_mat_names = [m[1:end-4] for m in list_o_mats]
    mats = [matrix_IO.load_mat(load_dir*mat) for mat in list_o_mats]
    return sort(Dict(zip(list_o_mat_names, mats)))
  end
end


## Examples
### Literature Examples
function generate_lit_cp_mats()
    lit_cp_mats = Dict()
    #21BS-M: Berman and Shaked-Monderer. “Copositive and Completely Positive Matrices.” (2021): p265,p302
    lit_cp_mats["ex01_n5_zd0.5_r5_ucpr0_21BS-Mp265"] = [1 1 0 0 1
                          1 2 1 0 0
                          0 1 3 1 0
                          0 0 1 4 1
                          1 0 0 1 5]

    lit_cp_mats["ex02_n5_zd0.5_r5_ucpr0_21BS-Mp302"] = [2 1 1 0 0
                      1 2 0 0 1
                      1 0 2 1 0
                      0 0 1 2 1
                      0 1 0 1 2]

    #98B: Barioli. “Completely positive matrices with a book-graph.”  Linear Algebra and its Applications 277 (1998): 11-31  
    lit_cp_mats["ex03_n5_zd0.5_r5_ucpr5_98Bex4.3"] =  [3 2  0 0 1
                                        2 5  6 0 0
                                        0 6 14 4 0
                                        0 0  4 9 1
                                        1 0  0 1 2]

    lit_cp_mats["ex04_n6_zd0.27_r6_ucpr0_98Bex1"] =[12  2 2  6 4 2
                          2 23 7  9 6 4
                          2  7 5  7 0 0
                          6  9 7 13 0 0
                          4  6 0  0 8 6
                          2  4 0  0 6 5]

    lit_cp_mats["ex05_n10_zd0.56_r10_ucpr0_98Bex2"] =[13  0 1 2 1 1 4 1 1 1
                                                      0 19 3 3 1 3 1 1 3 5
                                                      1  3 5 4 0 0 0 0 0 0
                                                      2  3 4 5 0 0 0 0 0 0
                                                      1  1 0 0 2 3 0 0 0 0
                                                      1  3 0 0 3 5 0 0 0 0
                                                      4  1 0 0 0 0 5 2 0 0
                                                      1  1 0 0 0 0 2 1 0 0
                                                      1  3 0 0 0 0 0 0 2 4
                                                      1  5 0 0 0 0 0 0 4 10]
      
    #11H: Haufmann. “On Completely Positive Matrices.” Master’s Thesis for the degree of Master of Applied Mathematics and Mechanics (2011): 30-31
    lit_cp_mats["ex06_n6_zd0.6_r6_ucpr0_11Hp30"] = [2 0 1  5  0  0
                                          0 1 1  0  0  0
                                          1 1 2  0  0  0
                                          5 0 0 31  5  6
                                          0 0 0  5 13  8
                                          0 0 0  6  8 17]

    #93SZ Salce, Zanardo “Completely Positive Matrices and Positivity of Least Squares Solutions.” (1993)  
    lit_cp_mats["ex07_n5_zd0.5_r5_ucpr0_93SZex2.3"] =  [ 2 1 0 0 1
                                              1 2 1 0 0
                                              0 1 2 1 0
                                              0 0 1 2 1
                                              1 0 0 1 3]

    #98XX  Xiang Xiang. “Notes on Completely Positive Matrices.” Linear Algebra and its Applications 271 (1998) 273-282  
    lit_cp_mats["ex08_n5_zd0.4_r4_ucpr6_98XXex4.1"] = [2 0 0 1 1
                      0 2 0 1 1
                      0 0 2 1 1
                      1 1 1 3 0
                      1 1 1 0 3]

    lit_cp_mats["ex09_n5_zd0.2_r5_ucpr9_98XXex4.2"] =[15   3  2 1 0
                                            3  18  7 0 2
                                            2   7 20 6 1
                                            1   0  6 9 2
                                            0   2  1 2 7]

    lit_cp_mats["ex10_n5_zd0.5_r5_ucpr7_98XXex4.3"] =[2 0 1 1 0
                          0 2 1 1 0
                          1 1 4 1 0
                          1 1 1 5 0
                          0 0 0 0 3]  

    #17BS-M:  Berman & Shaked-Monderer. “Completely Positive Matrices: Real, Rational, and Integral.” (2017) 
    lit_cp_mats["ex11_n5_zd0.0_r5_ucpr10_17BS-Mex1.6"] =[8 5 1 1 5
                                        5 8 5 1 1
                                        1 5 8 5 1
                                        1 1 5 8 5
                                        5 1 1 5 8]
      
    #10AD  Anstreicher,Dong. “A note on 5 x 5 Completely Positive Matrices” Linear Algebra and its Applications 433 (2010) 1001-1004 
    lit_cp_mats["ex12_n5_zd0.6_r4_ucpr0_10ADp4"] = [2 0 0 1 0
                                                    0 2 0 0 1
                                                    0 0 4 1 1
                                                    1 0 1 1 0
                                                    0 1 1 0 1]
    
    #04BX: Berman & Xu. “5 x 5 Completely positive matrices with a book-graph.”  Linear Algebra and its Applications 393 (2004): 55-71
    lit_cp_mats["ex13_n5_zd0.1_r3_ucpr0_04BXex4.1"] = [6 5 5 1 2
                                                      5 6 5 2 1
                                                      5 5 9 2 2
                                                      1 2 2 1 0
                                                      2 1 2 0 1]

    #16CS: Cedolin & Salce. “Completely positive matrices of order 5 with CP-graph.” Note Mat. 36 (2016) no. 1, 123-132
    lit_cp_mats["ex14_n5_zd0.4_r5_ucpr0_16CSex3.2"] =  [10  5  2  0  5
                                          5 10  5  0  0
                                          2  5 10  5  0
                                          0  0  5 10  5
                                          5  0  0  5 10]
   
    lit_cp_mats["ex15_n5_zd0.3_r5_ucpr0_16CSex3.3"] =[10  5   2  0    5
                                            5 10   5 10/3  0
                                            2  5  10  5    0
                                            0 10/3 5 10    5
                                            5  0   0  5   10]
        
    #21BN Bot & Nguyen. “Factorization of completely positive matrices using iterative projected gradient steps.” (2021) and P. Groetzner, M. Dur, A factorization method for completely positive matrices, Linear Algebra Appl. (2020).
    # So W, Xu C. A simple sufficient condition for complete positivity. Operat Matrices. 2015;9(1):233–9.
    lit_cp_mats["ex16_n5_zd0.0_r3_ucpr3_21BNeq81"] =[41 43  80  56 50
                                                    43 62  89  78 51
                                                    80 89 162 120 93
                                                    56 78 120 104 62
                                                    50 51  93  62 65] 

     #09BSA Burer Antreicher Duer. “The difference between 5 x 5 doubly nonnegative and completely positive matrices.” Linear Algebra and its Applications 431(9) (2009) 1539-1552
    lit_cp_mats["ex17_n5_zd0.5_r5_ucpr0_09BADp5"] = [24 11   0   0 11
                                                     11 24  11   0  0
                                                      0 11 124  11  0
                                                      0  0  11  24 11
                                                     11  0   0  11 24]

    #21PS Pfeffer & Samper. “The cone of 5 x 5 Completely Positive Matrices.”  (2021)
    lit_cp_mats["ex18_n5_zd0.0_r5_ucpr5_21PSp12"] = [53 32  1  4  26
                                                    32 30  5  2   5
                                                    1  5 11 10   3
                                                    4  2 10 21  17
                                                    26  5  3 17  35]

    #14N 2014, Nie, The A-Truncated K-Moment Problem
    lit_cp_mats["ex19_n5_zd0.1_r5_ucpr5_14N1"] =[6 4 1 2 2
                                                 4 6 0 1 3
                                                 1 0 3 1 2
                                                 2 1 1 2 1
                                                 2 3 2 1 5] 
    #03BS-M from Completely positive matrices, Berman Shaked-Monderer,  World scientific
    # lit_cp_mats["ex20_n6_zd2.0_r2_ucpr0_03BS-Mexe2.1"] = [  8  12 16 4 6  8
    #                                                         12 20 28 6 10 14
    #                                                         16 28 40 8 14 20
    #                                                         4   6  8 2  3  4
    #                                                         6  10 14 3  5  7
    #                                                         8  14 20 4  7 10]
    # Groetzner. “A method for Completely Positive and Nonnegative Matrix Factorization.” (2018)
    # lit_cp_mats["ex21_n7_zd2.0_r7_ucpr10_18Gp104"] = [13.7162 8.1090 10.0752  7.3940 10.7332 4.2551 9.8380
    #                                                   8.1090 9.9804  6.7089  6.1420  8.3044 3.9288 6.7772
    #                                                 10.0752 6.7089 69.6070  5.6133  8.3782 3.6309 6.6191
    #                                                   7.3940 6.1420  5.6133 10.2718  8.1573 4.3024 5.9480
    #                                                 10.7332 8.3044  8.3782  8.1573 12.8697 4.1992 9.1688
    #                                                   4.2551 3.9286  3.6309  4.3024  4.1092 2.5809 3.3640
    #                                                   9.8380 6.7772  6.6191  5.9480  9.1688 3.3640 13.1836]
    
    
    #14BSU: Matrices from:  I.M. Bomze, W. Schachinger, R. Ullrich. From seven to eleven: Completely positive matrices with high cp-rank. Linear Algebra and its Applications 459 (2014), 208 – 221.
    # lit_cp_mats["ex22_n7_zd2.0_r7_ucpr14_14BSUex1"] =  circmatrix(([531.0, 81, 150, 150, 81, 531, 926]))
    # lit_cp_mats["ex23_n7_zd2.0_r7_ucpr14#_14BSUtm7"] =  circmatrix(([108, 27, 4, 4, 27, 108, 163.0]))
    lit_cp_mats["ex20_n8_zd0.036_r8_ucpr18#_14BSUmex3"] =  [541.0 880 363 24 55 11 24 0
                                                          880 2007 1496 363 48 22 22 24
                                                          363 1496 2223 1452 363 24 22 11
                                                          24 363 1452 2325 1584 363 48 55
                                                          55 48 363 1584 2325 1452 363 24
                                                          11 22 24 363 1452 2223 1496 363
                                                          24 22 22 48 363 1496 2007 880
                                                          0 24 11 55 24 363 880 541]

    # lit_cp_mats["ex25_n9_zd2.0_r9_ucpr26#_14BSUmex2"] =  [2548 1628 363 60 55 55 60 363 1628
                                                          # 1628 2548 1628 363 60 55 55 60 363
                                                          # 363 1628 2483 1562 363 42 22 55 60
                                                          # 60 363 1562 2476 1628 363 42 55 55
                                                          # 55 60 363 1628 2548 1628 363 60 55
                                                          # 55 55 42 363 1628 2476 1562 363 60
                                                          # 60 55 22 42 363 1562 2483 1628 363
                                                          # 363 60 55 55 60 363 1628 2548 1628
                                                          # 1628 363 60 55 55 60 363 1628 2548]

    lit_cp_mats["ex21_n11_zd0.2_r11_ucpr32#_14BSUex5"] =  [781 0 72 36 228 320 240 228 36 96 0;
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
function generate_lit_cp_mats(save_dir::String)
  lit_cp_mats = generate_lit_cp_mats()
  for k in keys(lit_cp_mats)
    mat = lit_cp_mats[k]
    matrix_IO.save_mat(mat, save_dir*k*"_.csv")
  end
end

### Literature non Examples
function generate_lit_non_cp_mats()
  lit_non_cp_mats = Dict()
  ## 2014, Nie, The A-Truncated K-Moment Problem
  lit_non_cp_mats[1] = [1 1 0 0 1
                        1 2 1 0 0
                        0 1 2 1 0
                        0 0 1 2 1
                        1 0 0 1 6] 
  # Barioli. “Completely positive matrices with a book-graph.”  Linear Algebra and its Applications 277 (1998): 11-31 
  lit_non_cp_mats[2] = [7  1  2  2  1  1
                        1 12  1  3  3  5
                        2  1  2  3  0  0
                        2  3  3  5  0  0
                        1  3  0  0  2  4
                        1  5  0  0  4 10]
    #93SZ Salce, Zanardo “Completely Positive Matrices and Positivity of Least Squares Solutions.” (1993)  
    lit_non_cp_mats["n5_r5_ucpr0_93SZex2.3"] =  [ 1 1 0 0 1
                                                  1 2 1 0 0
                                                  0 1 2 1 0
                                                  0 0 1 2 1
                                                  1 0 0 1 3]                  



  for k in keys(lit_cp_mats)
    lit_cp_mats[k] = makediagone(lit_cp_mats[k])
  end                      
  return lit_cp_mats
end
function generate_lit_non_cp_mats(save_dir::String)
  lit_cp_mats = get_lit_non_cp_mats()
  for k in keys(lit_cp_mats)
    mat = lit_cp_mats[k]
    matrix_IO.save_mat(mat, save_dir*k*"_.csv")
  end
end


### Randomly generated Examples
# Generates a batch of cp-matrices, stores them, their construciton, and rank 
function generate_random_cp_mats(n_range=[5:9...], p_seq=[0.25:0.075:0.95...], save_dir=proj_dir)
  Random.seed!(373)
  ex = 1
  for n in n_range
      for p in p_seq 
          ext = lpad(ex,2,'0')
          ucpr, fact, M = get_random_cp_mat(n,p)
          nzd = rpad(round(get_nz_dense(M),digits=3),4,'0')
          r = LinearAlgebra.rank(M .+ 0.0)
          M_name = "ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
          matrix_IO.save_mat(M, save_dir*M_name)
          fact_name = "R_ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
          matrix_IO.save_mat(fact, save_dir*fact_name)
          ex = ex + 1
      end
  end
end
get_random_cp_mat(n::Int,p::Float64=0.5) = get_random_cp_mat(gen_random_sparcity_mat(n,p))
function get_random_cp_mat(M::Matrix{Int64})
    n = size(M)[1]
    nze = get_nonzero_entries(M)
    R = cat([ reᵢⱼ(n,e[1],e[2]) for e in nze]...,dims=2)
    rp = size(R)[2]
    return rp, R, makediagone(R*R' .+ 0.0)
end
gen_random_sparcity_mat(n,p=0.5) = Int.(Symmetric(map(x -> bernoulli(p), ones(n,n))) + I(n) .> 0)
reᵢⱼ(n::Int,i::Int,j::Int) = [k ∈ [i,j] ? rand() : 0 for k in 1:n]
get_nonzero_entries(M) = [(i,j) for i in 1:(size(M)[1]) for j in i:size(M)[1] if M[i,j] != 0]
bernoulli(p) =  rand() <= p ? 1 : 0
get_nz_dense(M) = round(1 - sum(M .== 0) / ((size(M)[1])*(size(M)[1]-1)),digits=2)

function get_random_cp_mat(n::Int,r::Int)
  a = rand(n,r)
  return r, a, makediagone(sum([a[:,k]*a[:,k]' for k ∈ 1:r ]))   
end


#Randomly generated Examples with bipartite graph support---------------------------------------------------------------------------
function gen_bipart_supp_mats(n_range, save_dir::String)
    ex = 0
    for n ∈ n_range
        m, a = gen_bipart_supp_mats(n)    
        ex += 1
        ext = lpad(ex, 2, '0')
        nzd = get_nz_dense(m)
        r = rank(m)
        ucpr = length(a)

        M_name = "ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
        matrix_IO.save_mat(m, save_dir*M_name)

        fact_name = "R_ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
        matrix_IO.save_mat(a, save_dir*fact_name)
    end

end

function gen_bipart_supp_mats(n) 
  M_supp = bipar_supp_mat(n)
  k = iseven(n) ? n^2 : n^2 + n
  return flesh_sup_mat(M_supp,k)
end
bipar_supp_mat(n) = [ bipar_check(i,j) ? 1 : 0 for i ∈ 1:n, j ∈ 1:n]
bipar_check(i,j) = (sum(iseven.([i,j])), sum(isodd.([i,j]))) == (1,1) || (i == j)


#---------------------------------------------------------------------------
function gen_asc_nzd_cp_mats(n::Int64, k::Int64, save_dir::String)
    asc_nzd_cp_mats = gen_asc_nzd_cp_mats(n,k)
    count = 0
    for (m, a) ∈ asc_nzd_cp_mats
        count += 1
        ext = lpad(count, 2, '0')
        nzd = get_nz_dense(m)
        r = rank(m)
        ucpr = length(a)

        M_name = "ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
        matrix_IO.save_mat(m, save_dir*M_name)

        fact_name = "R_ex$(ext)_n$(n)_nzd$(nzd)_r$(r)_ucpr$(ucpr)_.csv"
        matrix_IO.save_mat(a, save_dir*fact_name)
    end
end
gen_asc_nzd_cp_mats(n::Int64,k::Int64) = [flesh_sup_mat(M,k) for M in gen_asc_nzd_supp_mats(n)]
gen_asc_nzd_cp_mats(n::Int64,k::Int64,l::Int64) = flesh_sup_mat(gen_asc_nzd_supp_mats(n,l),k)


function flesh_sup_mat(M,k=size(M)[1])
  n = size(M)[1]
  mc = moments.get_maximal_cliques(M) 
  p = length(mc)
  K = rand_break_up(k,p)
  a_s = [embed(round.(rand(length(mc[k])),digits=2),mc[k],n) for k in 1:p for _ ∈ 1:K[k] ]
  return sum([a*a' for a in a_s]) , a_s
end
embed(v,β,n) = [i ∈ β ? popfirst!(v) : 0.0  for i in 1:n]
function rand_break_up(k,p)
  nar = ones(Int,p)
  while sum(nar) < k
      nar = nar + shuffle(moments.eᵢ(p,1))
  end
  return nar
end

gen_asc_nzd_supp_mats(n) = [gen_asc_nzd_supp_mats(n,l) for l = 1:combi(n,2)]
gen_asc_nzd_supp_mats(n,l) = zv2mat(n,get_rand_0_1_vec(l, combi(n,2)))
combi(n,k) = div(factorial(n),(factorial(n-k)*factorial(k)))
get_rand_0_1_vec(n1,n) =  shuffle([ones(Int64, n1)...,zeros(Int64,n-n1)...])
function zv2mat(n,zv)
    M = ones(Int,n,n)
    for i in 1:n, j in (i+1):n
        M[i,j] =  M[j,i] = popfirst!(zv)
    end
    return M
end

#---------------------------------------------------------------------------
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
#---------------------------------------------------------------------------
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