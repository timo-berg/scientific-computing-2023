using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("plot_def.jl")

###### Works but spectral radius is always 1 (apart from numerical rounding errors)
###### Something in the B_CGC matrix is probably wrong 


# c_values = 0:2:20
# n_values = 100:10:300

n_values = [100, 200, 300]
c_values = [0, 10, 50, 100]
spectral_radii = zeros(length(c_values), length(n_values))


for (i, c) in enumerate(c_values)
    for (j, n) in enumerate(n_values)
        A = get_A(n, c)

        B_CGC = get_B_CGC(A)
        B_eigvals = eigvals(B_CGC)
        ρ = maximum(abs.(B_eigvals))

        spectral_radii[i, j] = ρ
    end
end
println("maximum spectral radius: ", maximum(spectral_radii), " Minimum: ", minimum(spectral_radii))

# Heatmap of condition numbers
n_indices = 1:length(n_values)
c_indices = 1:length(c_values)

heatmap(n_indices, c_indices, round.(spectral_radii, digits=10), title="\n" * L"$\rho(B_{CGC})$ for different values of N and c" * "\n", xlabel="N", ylabel="c", color=:viridis, xticks=(n_indices, n_values), yticks=(c_indices, c_values))

savefig("plots/task_9_spectral_radius.png")