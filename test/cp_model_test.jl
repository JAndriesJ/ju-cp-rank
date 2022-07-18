module cp_model_test

using Test
using JuMP
src_dir = dirname(dirname(@__FILE__))*"\\src\\"
include(src_dir*"moments.jl")
include(src_dir*"matrix_IO.jl")
include(src_dir*"cp_model.jl")
const cpmo = cp_model

assets_dir = dirname(dirname(src_dir))*"\\assets\\"
data_dir = assets_dir*"data\\lit\\"

@testset "dense constraint sizes" begin
    mats = matrix_IO.load_mats(data_dir)
    for key in keys(mats)
        A = mats[key]
        n = size(A)[1]
        t = rand(2:4)
        E = length(moments.get_nonzero_entries(A)) - n

        model = JuMP.Model()
        @variable(model, Y[moments.make_mon_expo(size(A)[1], 2t)])

        sec_ord   = cpmo.make_sec_ord_mom_con(A,Y)
        dag_con   = cpmo.make_dag_con(A,t,Y)
        ddag_con  = cpmo.make_dag_con(A,t,Y,true)
        loc_con   = cpmo.make_loc_con(A,t,Y)
        xx_con    = cpmo.make_xx_con(A,t,Y)
        G_con     = cpmo.make_G_con(A,t,Y)
        
        @test length(sec_ord)   ==  binomial(n+1,2)
        @test length(dag_con)   == E
        @test length(dag_con[1]) ≤ binomial(n+2t-2,2t-2)
        @test length(ddag_con)  == length(dag_con) + n + 1
        @test length(loc_con)   == E + n
        @test length(xx_con)    == E 
        @test size(G_con)[1]    ≤ binomial(n+t-1,t-1)*n 
    end
end

@testset "sparse constraint sizes" begin
    mats = matrix_IO.load_mats(data_dir)
    for key in keys(mats)
        # key = [keys(mats)...][1]
        M = mats[key]
        n = size(M)[1]
        t = rand(2:4)
        E = length(moments.get_nonzero_entries(M)) - n

        model = JuMP.Model()
        mc    = moments.get_maximal_cliques(M)
        p     = length(mc)

        spar_inds = cp_model.make_spar_inds(mc,t)
        @variable(model, Y[spar_inds])
        spar_mom_mat_con     = cpmo.make_spar_mom_mat_con(M,t,Y)
        spar_sec_ord_mom_con = cpmo.make_spar_sec_ord_mom_con(M,Y) 
        spar_dag_con         = cp_model.make_spar_dag_con(M,t,Y,false)
        spar_ddag_con        = cp_model.make_spar_dag_con(M,t,Y,true)
        spar_xx_con          = cp_model.make_spar_xx_con(M,t,Y)
        spar_loc_con         = cpmo.make_spar_loc_con(M,t,Y)
        spar_G_con           = cpmo.make_spar_G_con(M,t,Y)
        spar_weak_G_con      = cpmo.make_spar_G_con(M,t,Y,true)
        obj_Y                = cpmo.make_obj(mc,Y)
    
        sumVₖ       =  sum([length(mc[i]) for i ∈ 1:p])
        sumVₖchoose2 = sum([binomial(length(mc[i]),2) for i ∈ 1:p])

        @test length(spar_mom_mat_con) == p 
        @test length(spar_sec_ord_mom_con) == E+n
        @test length(spar_dag_con) == sumVₖchoose2
        @test length(spar_ddag_con) == p + sumVₖ + sumVₖchoose2
        @test length(spar_xx_con) == sumVₖchoose2 
        @test length(spar_loc_con) == sumVₖ  + sumVₖchoose2
        @test length(spar_G_con) == p
        @test length(spar_weak_G_con) == p

        for k ∈ 1:p
            Vₖ = length(mc[k])
            Eₖ = length(moments.get_edges(mc[k],false))
            @test size(spar_mom_mat_con[k])[1] == binomial(Vₖ+t,t)  
            @test size(spar_xx_con[k])[1] == binomial(Eₖ+t-1,t-1)
            @test size(spar_G_con[k])[1] == n*binomial(Vₖ+t-1,t-1)
            @test size(spar_weak_G_con[k])[1] == Vₖ*binomial(Vₖ+t-1,t-1)
        end
        for k ∈ 1:p
            for e ∈ moments.get_edges(sort(mc[k]))
                size(spar_loc_con[k])[1] == binomial(length(mc[k])+t-1,t-1)
            end
        end
    end
end

@testset "cp_model.get_ξₜᶜᵖ t =1 " begin
    M =    [1 0 0 0 1
            0 2 0 0 1
            0 0 3 1 0
            0 0 1 4 1
            1 1 0 1 5] .+ 0.0

    t = 1
    ξₜᶜᵖⁱᵈ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"");
    ξₜᶜᵖⁱᵈ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"G");

    @test objective_value(ξₜᶜᵖⁱᵈ1) == 3.5581388700291354
    @test objective_value(ξₜᶜᵖⁱᵈ2) == 3.5581389111007464
    
    ξₜᶜᵖʷˢᵖ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"");
    ξₜᶜᵖʷˢᵖ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"G");
    
    @test objective_value(ξₜᶜᵖʷˢᵖ1) == 4.000000009870513
    @test objective_value(ξₜᶜᵖʷˢᵖ2) == 4.000000006724337

    ξₜᶜᵖˢᵖ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"");
    ξₜᶜᵖˢᵖ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"G");

    objective_value(ξₜᶜᵖˢᵖ1) == 4.000000009870513
    objective_value(ξₜᶜᵖˢᵖ2) == 4.000000027034827
end

@testset "cp_model.get_ξₜᶜᵖ t =2 " begin
    M =    [1 0 0 0 1
            0 2 0 0 1
            0 0 3 1 0
            0 0 1 4 1
            1 1 0 1 5] .+ 0.0
    t = 2
    ξₜᶜᵖⁱᵈ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"");
    ξₜᶜᵖⁱᵈ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"G");
    ξₜᶜᵖⁱᵈ3 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gdag");
    ξₜᶜᵖⁱᵈ4 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gddag");
    ξₜᶜᵖⁱᵈ5 , _ = cp_model.get_ξₜᶜᵖ(M,t,"id"*"Gddagxx");

    @test objective_value(ξₜᶜᵖⁱᵈ1) == 4.284770679528348
    @test objective_value(ξₜᶜᵖⁱᵈ2) == 4.994713677680638
    @test objective_value(ξₜᶜᵖⁱᵈ3) == 4.9971313589576445
    @test objective_value(ξₜᶜᵖⁱᵈ4) == 4.99843390279713
    @test objective_value(ξₜᶜᵖⁱᵈ5) == 4.99844411015471

    ξₜᶜᵖʷˢᵖ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"");
    ξₜᶜᵖʷˢᵖ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"G");
    ξₜᶜᵖʷˢᵖ3 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gdag");
    ξₜᶜᵖʷˢᵖ4 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gddag");
    ξₜᶜᵖʷˢᵖ5 , _ = cp_model.get_ξₜᶜᵖ(M,t,"wsp"*"Gddagxx");

    @test objective_value(ξₜᶜᵖʷˢᵖ1) == 4.563978851325215
    @test objective_value(ξₜᶜᵖʷˢᵖ2) == 4.596491388691128
    @test objective_value(ξₜᶜᵖʷˢᵖ3) == 4.596491857319117
    @test objective_value(ξₜᶜᵖʷˢᵖ4) == 4.596491420561562
    @test objective_value(ξₜᶜᵖʷˢᵖ5) == 4.596491350456679

    ξₜᶜᵖˢᵖ1 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"");
    ξₜᶜᵖˢᵖ2 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"G");
    ξₜᶜᵖˢᵖ3 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gdag");
    ξₜᶜᵖˢᵖ4 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gddag");
    ξₜᶜᵖˢᵖ5 , _ = cp_model.get_ξₜᶜᵖ(M,t,"sp"*"Gddagxx");   

    @test objective_value(ξₜᶜᵖˢᵖ1) == 4.563978851325215
    @test objective_value(ξₜᶜᵖˢᵖ2) == 5.000000063914473
    @test objective_value(ξₜᶜᵖˢᵖ3) == 5.000000016302208
    @test objective_value(ξₜᶜᵖˢᵖ4) == 5.000000071850989
    @test objective_value(ξₜᶜᵖˢᵖ5) == 5.000000069380393
end

end