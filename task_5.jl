# Compute and tabulate the condition number of the Symmetric Gauss-Seidel precondi-
# tioned operator M^{−1}_{SGS}A.

# TODO: Check if the large condition numbers for c = 10 are correct
# M_SGS = (D-E) * D^{-1} * (D - F)

using LinearAlgebra
using Plots
using LaTeXStrings
include("gauss_seidel.jl")
include("utils.jl")



function compare_condition_numbers(get_inv_M, method)
    c_values = 0:2:20
    n_values = 10:10:200
    cond_values = zeros(length(c_values), length(n_values))

    for (i, c) in enumerate(c_values)
        for (j, n) in enumerate(n_values)
            A = get_A(n, c)
            M_inv_SGS = get_inv_M(A)
            cond_nr = cond(M_inv_SGS * A)
            cond_values[i, j] = cond_nr
        end
    end

    # Heatmap of condition numbers
    p = heatmap(n_values, c_values, cond_values, title=L"Condition number of $M_{%$method}^{-1}A$", xlabel="N", ylabel="c", color=:viridis)
    return p
end


compare_condition_numbers(get_inv_M_SGS, "SGS")