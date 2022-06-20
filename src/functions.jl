using ArnoldiMethod

function energyspectrum(H, n::Int64, ϵ)
    decomp, history = partialschur(H, nev=n, which=SR(), tol=ϵ)

    E = decomp.eigenvalues[n] |> real

    E
end

function energyspectrum_BIG(H, n::Int64, ϵ)
    decomp, history = partialschur(H, nev=n, which=SR(), tol=ϵ)

    E = decomp.eigenvalues[n] |> real

    return BigFloat.(E)
end


"""
Function for sorting string in 1,2,3,4,5,6,7,8,9,10,11 ... 
"""
function natural(x, y)
    k(x) = [occursin(r"\d+", s) ? parse(Int, s) : s
            for s in split(replace(x, r"\d+" => s -> " $s "))]
    A = k(x)
    B = k(y)
    for (a, b) in zip(A, B)
        if !isequal(a, b)
            return typeof(a) <: typeof(b) ? isless(a, b) :
                   isa(a, Int) ? true : false
        end
    end
    return length(A) < length(B)
end