module matrix_IO_test
using Test

src_path = dirname(dirname(@__FILE__))*"\\src\\"
include(src_path*"matrix_IO.jl")
include(src_path*"cp_model.jl")

M = [1 0 0 0 1
     0 2 0 0 1
     0 0 3 1 0
     0 0 1 4 1
     1 1 0 1 5] .+ 0.0
     
_, ext_mom_dens = cp_model.get_ξₜᶜᵖ(M,1,"id"*""); 
_, ext_mom_spars = cp_model.get_ξₜᶜᵖ(M,2,"sp"*"Gdagxx"); 


test_dir = dirname(@__FILE__)*"\\"
@testset "save & load matrices" begin
    M = reshape(1:30,5,6)
    save_path = load_path = test_dir*"test_mat.csv"
    matrix_IO.save_mat(M, save_path)
    @test M == matrix_IO.load_mat(load_path) 

    src_dir = dirname(@__FILE__)*"\\"
    assets_dir = dirname(dirname(src_dir))*"\\assets\\"
    data_dir = assets_dir*"data\\literature\\"
    data_matrix_dict = matrix_IO.load_mats(data_dir)
    dmdk = data_matrix_dict.keys
    dmdv = data_matrix_dict.vals

    @test @assert dmdk[3] == "ex03_n5_nzd0.5_r5_ucpr5_98Bex4.3_"
    @test @assert dmdk[9] =="ex09_n5_nzd0.4_r4_ucpr0_10ADp4_"
    @test @assert dmdk[14] =="ex14_n5_nzd0.5_r5_ucpr0_09BADp5_"
    @test @assert dmdk[20] =="ex20_n10_nzd0.44_r10_ucpr0_98Bex2_"


    @test @assert isapprox(dmdv[12], [1.0  0.5       0.2  0.0       0.5
                                        0.5  1.0       0.5  0.333333  0.0
                                        0.2  0.5       1.0  0.5       0.0
                                        0.0  0.333333  0.5  1.0       0.5
                                        0.5  0.0       0.0  0.5       1.0],atol=0.000001)

    @test @assert isapprox(dmdv[17], [1.0       0.120386  0.258199  0.480384  0.408248  0.258199
                                        0.120386  1.0       0.652753  0.520483  0.442326  0.373002
                                        0.258199  0.652753  1.0       0.868243  0.0       0.0
                                        0.480384  0.520483  0.868243  1.0       0.0       0.0
                                        0.408248  0.442326  0.0       0.0       1.0       0.948683
                                        0.258199  0.373002  0.0       0.0       0.948683  1.0],atol=0.00001)

    @test @assert isapprox(dmdv[8], [ 1.0    0.625  0.125  0.125  0.625
                                        0.625  1.0    0.625  0.125  0.125
                                        0.125  0.625  1.0    0.625  0.125
                                        0.125  0.125  0.625  1.0    0.625
                                        0.625  0.125  0.125  0.625  1.0],atol=0.00001)   
                                        
    meta_data_frame = matrix_IO.get_mat_data(data_dir)
    @test names(meta_data_frame) == ["name", "ex", "n", "nzd", "r", "ucpr", "nc", "mc"]
    @test size(meta_data_frame) == (22, 8)
end

@testset "save & load moments" begin
    save_path = load_path = test_dir*"test_mom_dens.csv"
    matrix_IO.save_moments(ext_mom_dens, save_path) 
    @test Matrix(matrix_IO.load_moments(load_path)) == ext_mom_dens

    save_path = load_path = test_dir*"test_mom_sparse.csv"
    matrix_IO.save_moments(ext_mom_spars, save_path) 
    @test Matrix(matrix_IO.load_moments(load_path)) == ext_mom_spars
end



@testset "split_name" begin
    K = ["ex01_n5_nzd0.5_r5_ucpr0_21BS-Mp265_"
    "ex07_n5_nzd0.5_r5_ucpr7_98XXex4.3_"
    "ex14_n5_nzd0.5_r5_ucpr0_09BADp5_"
    "ex19_n8_nzd0.97_r8_ucpr18_14BSUmex3_"
    "ex22_n12_nzd0.73_r10_ucpr37_15BSUappc_"]
    @assert all((matrix_IO.split_name(K[1]) .== (5, 0.5, 5, 0, "01")))
    @assert all((matrix_IO.split_name(K[2]) .== (5, 0.5, 5, 7, "07")))
    @assert all((matrix_IO.split_name(K[3]) .== (5, 0.5, 5, 0, "14")))
    @assert all((matrix_IO.split_name(K[4]) .== (8, 0.97, 8, 18, "19")))
    @assert all((matrix_IO.split_name(K[5]) .== (12, 0.73, 10, 37, "22")))
end


@testset "load_mats" begin

end

end







