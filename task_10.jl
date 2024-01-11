using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")
include("plot_def.jl")

# TODO: spectral radius of what?
# TODO: solving the system directly after one step is a bit weird?

# Get problem
p = plot(title="Error of the Coarse Grid Correction Method", xlabel="Iteration", ylabel="Relative Error")
# ylims!(0, 1e-2)

c_values = [0]#[-5, 0, 5, 10, 100]
N_values = [10, 50, 100, 200, 300]
# c_values = [100]
for N in N_values
    for c in c_values
        A = get_A(N, c)
        b = get_b(N, construct_F(c))

        M_CGC_inv = get_M_CGC_inv_dirichilet(A)

        u, iters, residuals = iterative_solve(A, b, zeros(N - 1), M_CGC_inv, 1e-6, 1000, true)

        errors = [norm(r[1:end]) / norm(b) for r in residuals[1:end]]
        # Color them according to N value
        color = cgrad(:viridis)[Int(round((N - 0) / (maximum(N_values) - 0) * (length(cgrad(:viridis)) - 1)) + 1)]
        plot!(p, errors, label="N = $N", color=color)
    end
end

p



savefig("plots/task_10_cgc_convergence.png")