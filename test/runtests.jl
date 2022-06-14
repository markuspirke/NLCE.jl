using Test
using Revise
using NLCE

filename = "/Users/markuspirke/Desktop/Studium/Bachelorarbeit/Data_triangular/DreieckDTO16/BondFile5-4-em60-aut24.txt"
cluster = create_cluster(filename)
cluster_BIG = NLCE.create_cluster_BIG(filename)

typeof(cluster)
typeof(cluster_BIG)

tfim(cluster_BIG, 0.1, 1)

@testset "NLCE.jl" begin
    @test cluster.order == 6
    @test cluster.sites == 5
end