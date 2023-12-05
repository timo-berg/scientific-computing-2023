# Compute and tabulate the condition number of the Symmetric Gauss-Seidel precondi-
# tioned operator M^{−1}_{SGS}A.

# M_SGS = (D-E) * D^{-1} * (D - F)

using LinearAlgebra
using Plots
using LaTeXStrings
include("utils.jl")

function get_M_SGS(A)
    n = size(A)[1]
    D = Diagonal(A)
    E = -UpperTriangular(A)
    F = -LowerTriangular(A)
    return (D - E) * inv(D) * (D - F)
end


c_values = -100:10:20
n_values = 10:5:60
cond_values = zeros(length(c_values), length(n_values))

for c = c_values
    for n = n_values
        A = get_A(n, c)
        M_SGS = get_M_SGS(A)
        try
            cond_values[c_values.==c, n_values.==n] = cond(inv(M_SGS) * A)
        catch
            println("Error for c = $c, n = $n")
        end
    end
end

# Heatmap of condition numbers
heatmap(n_values, c_values, cond_values, title="Condition number of M_SGS^{-1}A", xlabel="N", ylabel="c", color=:viridis, legend=:none)