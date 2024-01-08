using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")

function get_projection_matrix(N, A)
    n_h, n_H = get_fine_and_coarse_nr_node(N)
    M_inv = get_M_CGC_inv(n_H, n_h, c)
    return I - A * M_inv
end

function solve_projected_CG(A, b, N)
    P = get_projection_matrix(N, A)
    A_projected = P * A
    b_projected = P * b
    init_guess = zeros(N - 1)
    u, residuals, iters_used = conjugate_gradient(A_projected, b_projected, init_guess, 1e-6, 1000)
    return u, residuals, iters_used
end

function task11_residual_plot()
    # Get problem
    N = 300
    c_values = -50:10:50
    p = plot(xlabel="Iteration", ylabel="Residual Error", yscale=:log10, title="Convergence of Projected CG")
    yticks!(10.0 .^(-5:1:4))
    for c in c_values
        A = get_A(N, c)
        b = get_b(N, construct_F(c))
        _, errors, iters_used = solve_projected_CG(A, b, N)
        plot!(0:iters_used, errors, label="N = $N, c = $c")
    end
    savefig(p, "plots/task11_residual_plot.pdf")
    p
end


# Get the system
task11_residual_plot()

# # Plot error
# p_error = plot(errors, title="Convergence of Multi-grid", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c", yscale=:log10)

# # Plot both
# plot(p, p_error, layout=(2, 1), size=(800, 800))