using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")
include("plot_def.jl")

function get_projection_matrix(N, A)
    M_inv = get_M_CGC_inv(A)
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
    N_values = [10, 100, 300]
    c_values = [50, 100]
    p = plot(xlabel="Iteration", ylabel="Log Residual Error", yscale=:log10, title="Convergence of Projected CG")
    yticks!(10.0 .^ (-5:1:4))
    for N in N_values
        for c in c_values
            A = get_A(N, c)
            b = get_b(N, construct_F(c))
            _, errors, iters_used = solve_projected_CG(A, b, N)
            plot!(0:iters_used, errors, label="N = $N, c = $c")
        end
    end
    p
end

task11_residual_plot()

savefig("plots/task_11_residual.png")
