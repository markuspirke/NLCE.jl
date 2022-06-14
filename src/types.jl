abstract type Cluster end

struct Cluster64 <: Cluster
    name::String
    order::Int64
    sites::Int64 #number of sites
    bonds::Vector{Tuple{Int64,Int64}} #bonds between sites
    C::Float64 #embedding factor
end

struct ClusterBIG <: Cluster
    name::String
    order::Int64
    sites::Int64 #number of sites
    bonds::Vector{Tuple{Int64,Int64}} #bonds between sites
    C::BigFloat #embedding factor
end