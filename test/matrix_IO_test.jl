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
end

@testset "save & load moments" begin
    save_path = load_path = test_dir*"test_mom_dens.csv"
    matrix_IO.save_moments(ext_mom_dens, save_path) 
    @test Matrix(matrix_IO.load_moments(load_path)) == ext_mom_dens

    save_path = load_path = test_dir*"test_mom_sparse.csv"
    matrix_IO.save_moments(ext_mom_spars, save_path) 
    @test Matrix(matrix_IO.load_moments(load_path)) == ext_mom_spars
end


@testset "load_mats" begin

end

end





