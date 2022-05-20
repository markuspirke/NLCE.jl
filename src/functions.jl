using ArnoldiMethod

function energyspectrum(H, n::Int64)
    decomp, history = partialschur(H, nev=n, which=SR())

    E = decomp.eigenvalues[n] |> real

    E
end

