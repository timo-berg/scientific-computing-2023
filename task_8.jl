using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")

###### Works but not sure if its correct? 

N = 300 # Nr of mesh elements. 
c = 60
A = get_A(N, c)
# n_h = N - 1 # Nr of internal nodes in fine grid
# n_H = Int(N / 2 - 1) # Nr of exernal nodes in coarse grid

B_CGC = get_B_CGC(A)

B_eigvals = eigvals(B_CGC)

# Absoulte values of eigenvalues in a plot
p = plot(title="\n" * L"Absolute Eigenvalues of $B_{CGC}$" * "\n")
abs_eigvals = map(x -> abs(x), B_eigvals)
scatter!(p, abs_eigvals, label="", ms=2, markerstrokewidth=0)
p

# savefig("plots/task_8_eigenvalues.png")


# # Plot solution
N = 10
A = get_A(N, c)
b = get_b(N, construct_F(c))
u_exact(x) = exp(x) * (1 - x)


M_CGC_inv = get_M_CGC_inv_dirichilet(A)
u, iters, residuals = iterative_solve(A, b, zeros(N - 1), M_CGC_inv, 1e-6, 10, true)

x = range(1 / N, 1 - 1 / N, length=N - 1)
p = plot(x, u, label="Numerical Solution", title="Coarse Grid Correction Iterative Solution", xlabel="x", ylabel="u(x)")
plot!(x, u_exact.(x), label="Exact Solution")

savefig("plots/task_8_cgc_solution_coarse.png")

M_CGC_inv = get_M_CGC_inv(A)
u, iters, residuals = iterative_solve(A, b, zeros(N - 1), M_CGC_inv, 1e-6, 10, true)

x = range(1 / N, 1 - 1 / N, length=N - 1)
p = plot(x, u, label="Numerical Solution", title="Coarse Grid Correction Iterative Solution", xlabel="x", ylabel="u(x)")
plot!(x, u_exact.(x), label="Exact Solution")

savefig("plots/task_8_cgc_solution_dirichilet.png")