using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")

###### Works but spectral radius is always 1 (apart from numerical rounding errors)
###### Something in the B_CGC matrix is probably wrong 


c_values = 0:2:20
n_values = 10:10:200
spectral_radii = zeros(length(c_values), length(n_values))


for (i, c) in enumerate(c_values)
    for (j, n) in enumerate(n_values)
        h = N - 1
        H = Int(N / 2 - 1)

        B_CGC = get_B_CGC(H, h, c)
        B_eigvals = eigvals(B_CGC)
        ρ = maximum(abs.(B_eigvals))

        spectral_radii[i, j] = ρ
    end
end

# Heatmap of condition numbers
heatmap(n_values, c_values, spectral_radii, title=L"$\rho(B_{CGC})$ for different values of $N$ and $c$", xlabel="N", ylabel="c", color=:viridis)