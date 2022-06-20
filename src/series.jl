using TaylorSeries, Elliptic
#series expansion up to order 16 for the triangular lattice ferromagnet
series_ferro(h0) = -2598117.904386449 * h0^16 - 660091.9352334255 * h0^15 - 170007.9698921716 * h0^14 - 44479.57261067488 * h0^13 - 11852.64218066848 * h0^12 - 3227.548748096352 * h0^11 - 902.046897887663 * h0^10 - 260.2234497071275 * h0^9 - 78.12725830073408 * h0^8 - 24.66796875 * h0^7 - 8.35546875 * h0^6 - 3.09375 * h0^5 - 1.359375 * h0^4 - 0.75 * h0^3 - 0.75 * h0^2
series_ferro_taylor(N) = taylor_expand(series_ferro, 0, order=N)

#analytic solution via elliptic integral for the TFIM chain 
λ(J) = 1 / J
x(λ) = 4λ / (1 + λ)^2
E_exact(J) = -J * 2 / π * (1 + λ(J)) * Elliptic.E(x(λ(J)))


