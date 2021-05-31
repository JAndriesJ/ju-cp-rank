using LinearAlgebra
const la = LinearAlgebra
using Test
include(dirname(dirname(@__FILE__))*"\\src\\cpMatrices.jl")
using .cpMatrices

cp_mats = cpMatrices.get_cp_mats()

@testset "cpMatrices" begin
    # Test the Size
    @test [size(cp_mats[key]) for key in keys(cp_mats)] == [(6, 6), (7, 7), (7, 7), (8, 8), (11, 11), (9, 9)]
    # Test the Rank
    @test [la.rank(cp_mats[key]) for key in keys(cp_mats)] ==  [2, 7, 7, 8, 11, 9]
    # Test the PSDness
    @test all([minimum(la.eigvals(cp_mats[key])) for key in keys(cp_mats)] .> -10^(-10))
    # Test nonegativity
    all([ la.transpose(cp_mats[key]) ≈ cp_mats[key] for key in keys(cp_mats)])
end
