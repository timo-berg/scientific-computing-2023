# PROBLEM: Real Eigenvalues have completely different scale than the Ritz values

using LinearAlgebra
using Plots
using LaTeXStrings
using ColorSchemes

include("utils.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")

# Choose n = 300, c = 60 based on previous results

# Get problem
N = 300
c = 60
A = get_A(N, c)
b = get_b(N, construct_F(c))

# Get Preconditioner
M_inv = get_inv_M_SGS(A)

# Plot eigenvalues of M_inv
# 

## Ritz values
# Get the Ritz values
x0 = zeros(N - 1)
tol = 1e-6
maxiter = 1000
residuals = preconditioned_cg_residuals(A, b, x0, tol, maxiter, M_inv)

# Normalize the residuals
residuals = normalize.(residuals)

# Get the Ritz values
p = plot()
color_palette = ColorSchemes.seaborn_rocket_gradient
R = []
for (i, residual) in enumerate(residuals)
    # Concatenate the residuals as columns to a matrix
    if i == 1
        global R = residual
    else
        global R = hcat(R, residual)

        # Get the Ritz values every 10 iterations
        if i % 10 == 0
            T = R' * A * R
            color = color_palette.colors[Int(round((i - 0) / (length(residuals) - 0) * (length(color_palette.colors) - 1)) + 1)]
            plot!(p, eigvals(T), label="", ms=3, color=color)
        end
    end
end

# Normalize eig_M_inv to eig_T
plot!(p, eigvals(M_inv), label="Real Eigenvalues", ms=3, color=:blue)
p


