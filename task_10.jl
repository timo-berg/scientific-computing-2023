using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")

# TODO: spectral radius of what?
# TODO: solving the system directly after one step is a bit weird?

# Get problem
N = 100
c = 1

n_h, n_H = get_fine_and_coarse_nr_node(N)

# Get the system
A = get_A(N, c)
b = get_b(N, construct_F(c))

M_CGC_inv = get_M_CGC_inv(n_H,n_h,c)
#M_J_inv = Diagonal(1 ./ diag(A))  # This works, just to sanity check the iterative solvers

u, iters = iterative_solve(A, b, zeros(N - 1), M_CGC_inv, 1e-6, 100)
println(iters, size(u), u[1:2])
println(norm(A*u - b))


# Plot solution
u_exact(x) = exp(x) * (1 - x)

x = range(1/N, 1-1/N, length=N - 1)
p = plot(x, u, label="Numerical Solution", title="CGC iterative Solution", xlabel="x", ylabel="u(x)")
plot!(x, u_exact.(x), label="Exact Solution")
xlims!(p,0.9, 1)

# # Plot error
# p_error = plot(errors, title="Convergence of Multi-grid", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c", yscale=:log10)

# # Plot both
# plot(p, p_error, layout=(2, 1), size=(800, 800))