module TestMoments

using Test
include(dirname(dirname(@__FILE__))*"\\src\\Moments.jl")
using .Moments

@testset "TestMoments" begin
    for i in 1:5
        n = rand(6:11)
        t = rand(2:4)

        mon_expoT = Moments.make_mon_expo(n, t, true)
        mon_expoF = Moments.make_mon_expo(n, t, false)
        mon_expo_matT = Moments.make_mon_expo_mat(n, (t,t), true)
        mon_expo_matF = Moments.make_mon_expo_mat(n, (t,t), false)
        mom_expo_keys = Moments.make_mom_expo_keys(n, t)
## Check dimensions of the entries.
        @test size(mon_expoT[1])[1]  == n
        @test size(mon_expoF[1])[1]  == n
        @test size(mon_expo_matT[1])[1]  == n
        @test size(mon_expo_matT[1])[1]  == n
        @test size(mom_expo_keys[1])[1]  == n
## Check the sizes of the arrays.
        @test size(mon_expoT)[1] == binomial(n+t,t)
        @test size(mon_expoF)[1] == binomial(n+t,t) - binomial(n+t-1,t-1)
        @test size(mon_expo_matT)[1] == binomial(n+t,t)
        @test size(mon_expo_matF)[1] == binomial(n+t,t) - binomial(n+t-1,t-1)
        @test size(mom_expo_keys)[1] == binomial(n + 2*t, 2*t)
    end
end

end
