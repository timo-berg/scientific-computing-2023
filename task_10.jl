using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")
include("plot_def.jl")

# TODO: spectral radius of what?
# TODO: solving the system directly after one step is a bit weird?

# Get problem
N = 20
c = 1

# Get the system
A = get_A(N, c)
b = get_b(N, construct_F(c))

M_CGC_inv = get_M_CGC_inv(A)

u, iters = iterative_solve(A, b, zeros(N - 1), M_CGC_inv, 1e-6, 10)


# Plot solution
u_exact(x) = exp(x) * (1 - x)

x = range(1 / N, 1 - 1 / N, length=N - 1)
p = plot(x, u, label="Numerical Solution", title="Coarse Grid Correction Iterative Solution", xlabel="x", ylabel="u(x)")
plot!(x, u_exact.(x), label="Exact Solution")

savefig("plots/task_10_cgc_solution.png")