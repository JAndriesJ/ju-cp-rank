module moments_test
using Test

src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"cp_matrices.jl")
include(src_dir*"nn_matrices.jl")
using .moments
const mo = moments

@testset "standerd basis vector eᵢ" begin
    for n in 6:13
        for k in 1:n
            e = eᵢ(n,k)
            @test length(e) == n
            @test maximum(e) == 1
            @test minimum(e) == 0
            @test sum(e) == 1
        end
    end
end

@testset "make_mon_expo sizes le and not le" begin
    for i in 1:5
        n = rand(6:11)
        t = rand(1:4)
        tt = (rand(1:4),rand(2:4))
        mon_expo_t_le  = mo.make_mon_expo(n,t; isle = true)
        mon_expo_tt_le = mo.make_mon_expo(n,tt; isle = true) 

        @test size(mon_expo_t_le)[1] == binomial(n+t,t)
        @test size(mon_expo_tt_le) == (binomial(n+tt[1],tt[1]), binomial(n+tt[2],tt[2]) )
        @test size(mon_expo_t_le[rand(1:binomial(n+t,t))])[1]  == n
        @test size(mon_expo_tt_le[rand(1:binomial(n+tt[1],tt[1])+binomial(n+tt[2],tt[2]))])[1]  == n

        mon_expo_t  = mo.make_mon_expo(n,t; isle = false)
        mon_expo_tt = mo.make_mon_expo(n,tt; isle = false) 
        @test size(mon_expo_t)[1] == binomial(n+t-1,n-1)
        @test size(mon_expo_tt) == (binomial(n+tt[1]-1,n-1),binomial(n+tt[2]-1,n-1))
    end
end

@testset "make_mon_expo(mon) sizes le and not le" begin
    for i in 1:5
        n = rand(6:11)
        t = rand(1:4)
        tt = (rand(1:4),rand(2:4))

        mom = mo.make_mon_expo(n,t)
        mom₁= mo.make_mon_expo(n,tt[1])
        mom₂= mo.make_mon_expo(n,tt[2])

        mon_expo_mom =  mo.make_mon_expo(mom)
        mon_expo_mom_mom =  mo.make_mon_expo(mom₁, mom₂)

        @test mon_expo_mom == mo.make_mon_expo(n,(t,t))
        @test mon_expo_mom_mom == mo.make_mon_expo(n,tt)

        mom = mo.make_mon_expo(n,t; isle = true)
        mom₁= mo.make_mon_expo(n,tt[1]; isle = true)
        mom₂= mo.make_mon_expo(n,tt[2]; isle = true)

        mon_expo_mom =  mo.make_mon_expo(mom)
        mon_expo_mom_mom =  mo.make_mon_expo(mom₁, mom₂)

        @test mon_expo_mom == mo.make_mon_expo(n,(t,t); isle = true)
        @test mon_expo_mom_mom == mo.make_mon_expo(n,tt; isle = true)
    end
end
 
@testset "make_mon_expo misc" begin
    MonBase = mo.make_mon_expo(2,3,isle=false)
    @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]

    MonBase = mo.make_mon_expo(2,4,isle=false)
    @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]

    MB = mo.make_mon_expo(2,(3,3),isle=true)
    @test size(MB) == (binomial(2+3,3),binomial(2+3,3))

    @test mo.make_mon_expo(mo.make_mon_expo(2,3)) == mo.make_mon_expo(2,(3,3))
end

@testset "make_mon_expo sparse" begin
    # make_mon_expo(t::Int,M::Matrix{Float64},isnn=false) 
    # make_mon_expo(t::Int,M::Matrix{Float64},j) 
    # make_mon_expo(t::Int,M::Matrix{Float64},j::Int,isnn=false) 
end

@testset "get_monomial_cliques" begin
    for z in 1:3
        n = rand(6:11)
        t = rand(1:4)
        r = rand(n-1:binomial(n,2))
        M = cp_matrices.gen_random_supp_mat(n,r) .+ 0.0 
        mc = mo.get_maximal_cliques(M) 
        p = length(mc)
        mon_cliq_tt_M = mo.get_monomial_cliques((t,t), M)
        mon_cliq_t_M = mo.get_monomial_cliques(t,M)
        mon_cliq_t_mc = mo.get_monomial_cliques(n,t,mc) 
        @test length(mon_cliq_tt_M) == p + 1
        for i ∈ 1:p
            @test length(mon_cliq_t_mc[i]) == binomial(length(mc[i])+t,t)
            @test mo.get_monomial_cliques(t,M)[i] == mon_cliq_t_mc[i]
            @test mo.make_mon_expo(mon_cliq_t_mc[i],mon_cliq_t_mc[i]) == mon_cliq_tt_M[i]
        end
        for j = 1:n
            mcj = [c for c in mc if j ∈ Set(c)]
            mon_cliq_t_m_j = mo.get_monomial_cliques(t,M,j)
            nar = [ sum(m) for m in mon_cliq_t_m_j[1:end-1]]
            @test length(nar) == length(mcj)
            for i in length(nar)
                @test all(nar[i][mcj[i]] .> 0)
            end
        end
    end
end 

@testset "expo_kron M" begin
    for i ∈ 1:2
        n = rand(1:10)
        t = (rand(1:3),rand(1:4))
        A = mo.make_mon_expo(n,t)
        B = mo.make_mon_expo(n,t)
        AOXB = mo.expo_kron(A,B)
        @test size(AOXB) == size(A) .* size(B) 
        @test AOXB[1,1] == A[1,1] + B[1,1]
    end
end

@testset "Misc" begin
    M = [1 0 ; 0 1]
    @test mo.get_zero_entries(M) == [(1, 2)]
    @test mo.get_nonzero_entries(M) == [(1, 1), (2, 2)]

    M = [1 0 0 0 1
        0 2 0 0 1
        0 0 3 1 0
        0 0 1 4 1
        1 1 0 1 5] .+ 0.0

    mo.get_maximal_cliques(M) == [[4, 5], [2, 5], [1, 5], [3, 4]]
    mo.get_edges([1,2,3,4]) == [[1, 1],[1, 2],[1, 3],[1, 4],[2, 2],[2, 3],[2, 4],[3, 3],[3, 4],[4, 4]]

    n = rand(1:15)
    M = cp_matrices.gen_random_supp_mat(n,rand(1:binomial(n,2)))
    mc = mo.get_maximal_cliques(M) 
    for c in mc
        @test M[c,c] == ones(length(c),length(c))
    end
end

@testset "make_Gₘ_adj_mat" begin
    for i ∈ 1:5
        m,n = rand(1:150,2)
        M = nn_matrices.make_NN_mat(m,n)
        sqrtMₘₐₓ = sqrt(maximum(M))
        Gₘ_adj_mat = mo.make_Gₘ_adj_mat(M)
        @test size(Gₘ_adj_mat) == (m+n,m+n)
        @test unique(Gₘ_adj_mat) == [1.0, 0.0]
        @test Gₘ_adj_mat[1:m,1:m] == ones(m,m)
        @test Gₘ_adj_mat[m+1:m+n,m+1:m+n] == ones(n,n)
    end
end

end



            # Moments
            @testset "make_mon_expo(n::Int,t::Int, isLeq::Bool = true)" begin
                for n ∈ 3:9, t ∈ 0:4
                    @test size(make_mon_expo_arr(n,t))[2] == n
                    MonBase = make_mon_expo(n,t,true)
                    #@test length(MonBase) == che(n,t)
                end
                for n ∈ 3:9, t ∈ 0:3
                    MonBase = make_mon_expo(n,t,true)
                    @test length(MonBase) == binomial(n +t,t)
                end

                MonBase = make_mon_expo(2,3,false)
                @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]

                MonBase = make_mon_expo(2,4,false)
                @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]
            end

            @testset "get_std_base_vec(n::Int,k::Int)" begin
                n =7
                for k in 1:n
                    e = get_std_base_vec(n,k)
                    @test length(e) == n
                    @test maximum(e) == 1
                    @test minimum(e) == 0
                    @test sum(e) == 1
                end
             end

            @testset "get_mon_index(B,α)" begin
                MonBase =  make_mon_expo(2,4,false)
                α = [2  2]
                @test  get_mon_index(MonBase, α) == 3

                MonBase = make_mon_expo(2,3)
                α = [2.0 1.0]
                @test  get_mon_index(MonBase, α) == 8


                MonBase = make_mon_expo(5,4)
                α = [0  0 0 2 2]
                @test  get_mon_index(MonBase, α) == 124
            end

            @testset "make_mom__expo_mat_dict(n::Int,t::Int)" begin

                for n ∈ 6:10 , t ∈ 2:4
                    local MonExp = make_mom_expo_mat_dict(n,t)
                    local MonBase = make_mon_expo(n,t)
                    for key in keys(MonExp)
                        for val in MonExp[key]
                            @test   key == MonBase[val[1]] + MonBase[val[2]]
                        end
                    end
                end
            end
        end