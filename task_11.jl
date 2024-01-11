using LinearAlgebra
using Plots
using LaTeXStrings
using Printf
include("multi_grid.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")
include("plot_def.jl")

function get_projection_matrix(A)
    M_inv = get_M_CGC_inv(A)
    return I - A * M_inv
end

function get_effective_condition_number(Ah)
    eigenvalues = eigvals(Ah)
    min_nonzero_eigenvalue = minimum(eigenvalues[eigenvalues .> 1e-8])
    return maximum(eigenvalues) / min_nonzero_eigenvalue 
end

function task11_residual_plot()
    # Get problem
    N_values = [10, 100, 300]
    c_values = [50, 100]
    p = plot(xlabel="Iteration", ylabel="Log Residual Error", yscale=:log10, title="Convergence of Projected CG")
    xticks!(p, 1:13)
    for N in N_values
        for c in c_values
            A = get_A(N, c)
            b = get_b(N, construct_F(c))
            P = get_projection_matrix(A)
            A_projected = Symmetric(P * A) # From analytical inspection this matrix is symmetric
            b_projected = P * b
            _, errors, iters_used, _ = conjugate_gradient(A_projected, b_projected, zeros(N-1), 1e-9, 2000)
            eff_cond = get_effective_condition_number(A_projected)
            plot!(errors, label="N = $N, c = $c, kappa = $(@sprintf("%.4f", eff_cond))")
        end
    end
    p
end

task11_residual_plot()

# savefig("plots/task_11_residual.png")
