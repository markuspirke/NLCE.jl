using LinearAlgebra, SparseArrays

"""
Takes as an input two tensors and gives as an output the tensorprodukt or kroneckerprodukt between
x and y 
"""
⊗(x, y) = kron(x, y)

struct Cluster
    name::String
    order::Int64
    sites::Int64 #number of sites
    bonds::Vector{Tuple{Int64,Int64}} #bonds between sites
    C::Float64 #embedding factor
end

"""
Cluster type with number sites, bonds, embedding factor and string name.
J interaction strength
magnet: +1 ferromagnet -1 antiferromagnet
"""
function tfim(C::Cluster, J::Float64, magnet::Int64)
    Z = [1 0; 0 -1] |> sparse
    X = [0 1; 1 0] |> sparse
    I = [1 0; 0 1] |> sparse

    N = C.sites

    interaction = fill(I, N)

    field = fill(I, N)
    field[1] = Z

    H = spzeros(Float64, 2^N, 2^N)

    for (i, j) in C.bonds

        interaction[i] = X
        interaction[j] = X

        H -= magnet * J * foldl(⊗, interaction)

        #interaction = fill(I, N)
        interaction[i] = I
        interaction[j] = I

    end

    for i in 1:N
        H -= magnet * foldl(⊗, field)
        field = circshift(field, 1)
    end
    H += N * foldl(⊗, interaction)

    H
end

"""
N numbers of sites. J strength of interaction between neighbouring sites. 
Creates H = J Σ XX' +  Σ Z with open boundary conditions. H is represented as a sparse matrix.
"""
function tfim_chain(N::Int64, J::Float64)
    Z = [1 0; 0 -1] |> sparse
    X = [0 1; 1 0] |> sparse
    I = [1 0; 0 1] |> sparse

    first_term_ops = fill(I, N)
    first_term_ops[1] = X
    first_term_ops[2] = X

    second_term_ops = fill(I, N)
    second_term_ops[1] = Z

    H = spzeros(Int, 2^N, 2^N)

    for i in 1:N-1
        H -= J * foldl(⊗, first_term_ops)


        first_term_ops = circshift(first_term_ops, 1)

    end

    for i in 1:N
        H -= foldl(⊗, second_term_ops)
        second_term_ops = circshift(second_term_ops, 1)
    end

    H
end


"""
N numbers of sites. J strength of interaction between neighbouring sites. 
Creates H = J Σ XX' +  Σ Z with periodic boundary conditions. H is represented as a sparse matrix.
"""
function tfim_chain_periodic(N::Int64, J::Float64)
    Z = [1 0; 0 -1] |> sparse
    X = [0 1; 1 0] |> sparse
    I = [1 0; 0 1] |> sparse

    first_term_ops = fill(I, N)
    first_term_ops[1] = X
    first_term_ops[2] = X

    second_term_ops = fill(I, N)
    second_term_ops[1] = Z

    H = spzeros(Int, 2^N, 2^N)

    for i in 1:N
        H -= J * foldl(⊗, first_term_ops)


        first_term_ops = circshift(first_term_ops, 1)

    end

    for i in 1:N
        H -= foldl(⊗, second_term_ops)
        second_term_ops = circshift(second_term_ops, 1)
    end

    H
end