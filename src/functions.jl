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

