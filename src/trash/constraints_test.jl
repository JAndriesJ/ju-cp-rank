module constraints_test
using Test
using JuMP
proj_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(proj_dir*"cpMatrices.jl")
include(proj_dir*"Moments.jl")
include(proj_dir*"Constraints.jl")
using .cpMatrices
using .Moments
using .Constraints

cp_mats = cpMatrices.get_cp_mats()
for key in keys(cp_mats)
    A = cp_mats[key ]
    n = size(A)[1]
    t = rand(2:4)

    model = JuMP.Model()
    list_of_keys = Moments.make_mom_expo_keys(n, t) # Define variables in the moment matrix.
    @variable(model, Lx[list_of_keys] )

    dag_con   = Constraints.make_dag_con(A,t,Lx)
    loc_con   = Constraints.make_loc_con(A,t,Lx)
    xx_con    = Constraints.make_xx_con(A,t,Lx)
    weakG_con = Constraints.make_weakG_con(A,t,Lx)
    G_con     = Constraints.make_G_con(A,t,Lx)

    @testset "Constraint sizes" begin
        @test length(dag_con) == n*(n+1)/2 + 1
        @test length(loc_con) == n*(n+1)/2
        @test length(xx_con) == n*(n+1)/2 - n
        @test length(weakG_con) == t - 1
        @test size(G_con)[1]  == binomial(n+t-1,t-1)*n
    end
end

function run_tests()
    @testset "sparse" begin
        # @test var_inds = get_var_inds(mc,t) 
        # @test sort(unique([ b  for (a,b) ∈ var_inds])) ==  sort(mom.make_mon_expo(length(mc[k]),2*t))
        # @test length(spar_pos_cons) == length(cat(get_maximal_cliques(M)...,dims=1))
        # @test length(spar_loc_cons) == length(cat([nz(m) for m in get_maximal_cliques(M)]...,dims=1))
    end
end


end 

