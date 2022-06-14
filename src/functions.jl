using ArnoldiMethod

function energyspectrum(H, n::Int64)
    decomp, history = partialschur(H, nev=n, which=SR())

    E = decomp.eigenvalues[n] |> real

    E
end

function energyspectrum_BIG(H, n::Int64)
    decomp, history = partialschur(H, nev=n, which=SR())

    E = decomp.eigenvalues[n] |> real

    return BigFloat(E)
end

function energyspectrum_BIG(H, n::Int64, tol)
    decomp, history = partialschur(H, nev=n, which=SR(), tol)

    E = decomp.eigenvalues[n] |> real

    return BigFloat(E)
end