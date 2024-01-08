using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")

# TODO: spectral radius of what?
# TODO: solving the system directly after one step is a bit weird?

# Get problem
N = 1000
c = 20


# Get the system
A = get_A(N, c)
b = get_b(N, construct_F(c))

# Solve and measure time
@time u, errors = multigrid(A, b, 1e-6, 1000, conjugate_gradient_solver)

# Compare against direct solution
# @time u_direct = A \ b



# Plot solution
u_exact(x) = exp(x) * (1 - x)

x = range(0, 1, length=N - 1)
p = plot(x, u, label="Numerical Solution", title="Multi-grid Solution", xlabel="x", ylabel="u(x)")
plot!(x, u_exact.(x), label="Exact Solution")

# Plot error
p_error = plot(errors, title="Convergence of Multi-grid", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c", yscale=:log10)

# Plot both
plot(p, p_error, layout=(2, 1), size=(800, 800))