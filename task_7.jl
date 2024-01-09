# PROBLEM: Real Eigenvalues have completely different scale than the Ritz values

using LinearAlgebra
using Plots
using LaTeXStrings
using ColorSchemes

include("utils.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")

# TODO: Ritz values are too big in magnitude. Why?

# Choose n = 300, c = 60 based on previous results

function plot_ritz_values(get_M_inv, method)
    # Get problem
    N = 300
    c = 60
    A = get_A(N, c)
    b = get_b(N, construct_F(c))

    # Get Preconditioner
    M_inv = get_M_inv(A)

    ## Ritz values
    ritz_sample_rate = method == "Gauss-Seidel" ? 10 : 1
    # Get the Ritz values
    x0 = zeros(N - 1)
    tol = 1e-6
    maxiter = 1000

    x, convergence, maxiter, residuals = preconditioned_cg(A, b, x0, tol, maxiter, M_inv)
    # x, convergence, maxiter, residuals = conjugate_gradient(M_inv * A, M_inv * b, x0, tol, maxiter)

    # Normalize the residuals
    residuals = normalize.(residuals)

    # Get the Ritz values
    p = plot(title="Ritz Values for $method")
    color_palette = reverse(ColorSchemes.seaborn_rocket_gradient)
    R = []
    for (i, residual) in enumerate(residuals)
        # Concatenate the residuals as columns to a matrix
        if i == 1
            global R = residual
        else
            global R = hcat(R, residual)

            # Get the Ritz values every 10 iterations
            if i % ritz_sample_rate == 0
                T = R' * M_inv * A * R
                color = color_palette.colors[Int(round((i - 0) / (length(residuals) - 0) * (length(color_palette.colors) - 1)) + 1)]

                # Plot the eigenvalues of T such that they are evenly distributed on the x axis            
                # scatter!(p, 0:1/i:1,abs.(eigvals(T)), label="", ms=1, color=color)
                plot!(p, abs.(eigvals(T)), label="", ms=3, color=color)
            end
        end
    end

    # Plot the real eigenvalues of M_inv * A
    plot!(p, abs.(eigvals(M_inv * A)), label="", ms=3, color=:blue)

    return p
end

plot_ritz_values(get_inv_M_SGS, "Gauss-Seidel")
