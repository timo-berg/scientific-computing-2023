# Compute and tabulate the condition number of the Symmetric Gauss-Seidel precondi-
# tioned operator M^{−1}_{SGS}A.

# M_SGS = (D-E) * D^{-1} * (D - F)

using LinearAlgebra
using Plots
using LaTeXStrings
include("utils.jl")

c_values = -20:2:20
n_values = 10:10:200
cond_values = zeros(length(c_values), length(n_values))

for (i, c) in enumerate(c_values)
    for (j, n) in enumerate(n_values)
        A = get_A(n, c)
        M_inv_SGS = get_inv_M_SGS(A)
        cond_nr = cond(M_inv_SGS * A)
        cond_values[i, j] = log(cond_nr)
    end
end

# Heatmap of condition numbers
heatmap(n_values, c_values, cond_values, title=L"Log Condition number of $M_{SGS}^{-1}A$", xlabel="N", ylabel="c", color=:viridis)