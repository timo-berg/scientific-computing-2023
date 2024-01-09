using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")

###### Works but spectral radius is always 1 (apart from numerical rounding errors)
###### Something in the B_CGC matrix is probably wrong 


c_values = 0:2:20
n_values = 100:10:300
spectral_radii = zeros(length(c_values), length(n_values))


for (i, c) in enumerate(c_values)
    for (j, n) in enumerate(n_values)
        h = n - 1
        H = Int(n / 2 - 1)
        A = get_A(n, c)

        B_CGC = get_B_CGC(H, h, c, A)
        B_eigvals = eigvals(B_CGC)
        ρ = maximum(abs.(B_eigvals))

        spectral_radii[i, j] = ρ
    end
end
println("maximum spectral radius: ", maximum(spectral_radii), " Minimum: ", minimum(spectral_radii))

# Heatmap of condition numbers
heatmap(n_values, c_values, spectral_radii, title=L"$\rho(B_{CGC})$ for different values of $N$ and $c$", xlabel="N", ylabel="c", color=:viridis)