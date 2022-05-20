using Test
using Revise
using NLCE

filename = "/Users/markuspirke/Desktop/Studium/Bachelorarbeit/Markus/Output/DreieckFullO12/BondFile5-4-em60-aut24.txt"
cluster = create_cluster(filename)

@testset "NLCE.jl" begin
    @test cluster.order == 6
    @test cluster.sites == 5
end