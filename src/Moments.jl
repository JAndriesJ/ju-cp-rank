module moments
using Test

export  eᵢ,
        eᵢⱼ,
        make_mon_expo,
        run_tests

"""The standard basis vector eᵢ in dimension n"""
eᵢ(n::Int,i::Int) = [Int(j==i) for j in 1:n]
eᵢⱼ(n::Int,i::Int,j::Int) = [k ∈ [i,j] ? 1 : 0 for k in 1:n]

"""[x]≦ₜ := [xᵅ for all α ∈ ℕⁿₜ] or [x]₌ₜ  := [xᵅ for all α ∈ ℕⁿ≤ₜ]"""
function make_mon_expo(n::Int,t::Int; isle::Bool = true)
    @assert typeof(t) == Int64
    t == 0 ? (return [eᵢ(n,0)]) : 0
    tmp = make_mon_expo(n,t-1;isle=isle)
    M_vec = reshape([m + eᵢ(n,i) for i ∈ 1:n, m ∈ tmp],:,1)
    return unique(isle ? vcat(tmp,M_vec) : M_vec)
end

"""[x]≦ₜ[x]ᵀ≦ₜ or [x]₌ₜ[x]ᵀ₌ₜ"""
function make_mon_expo(n::Int,t::Tuple{Int,Int}; isle::Bool = true)
    M_vec1      = make_mon_expo(n,t[1]; isle=isle)
    M_vec2      = make_mon_expo(n,t[2]; isle=isle)
    return [mi+mj for mi in M_vec1, mj in M_vec2]
end
        
###
function run_tests()
    @testset "standerd basis vector" begin
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

    @testset "make_mon_expo" begin
        for n ∈ 3:9, t ∈ 0:4
            @test size.(make_mon_expo(n,t,isle=true))[1] == (n,)
            @test length(make_mon_expo(n,t,isle=true)) == binomial(n+t,t)
        end

        MonBase = make_mon_expo(2,3,isle=false)
        @test MonBase ==  [ [3, 0],[2, 1],[1, 2],[0, 3]]

        MonBase = make_mon_expo(2,4,isle=false)
        @test MonBase ==  [[4, 0],[3, 1],[2, 2],[1, 3],[0, 4]]

        MB = make_mon_expo(2,(3,3),isle=true)
        @test size(MB) == (binomial(2+3,3),binomial(2+3,3))
    end
end

end
