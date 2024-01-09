using LinearAlgebra
using Plots
include("utils.jl")
include("gauss_seidel.jl")
include("multi_grid.jl")
include("task_2.jl")
include("task_3.jl")
include("task_4.jl")
include("plot_def.jl")
# Verify that the multigrid method works for the problem
function test_TGM()
    # Get problem
    N = 300
    c = 10

    n_h, n_H = get_fine_and_coarse_nr_node(N)

    # Get the system
    A = get_A(N, c)
    b = get_b(N, construct_F(c))

    M_TGM_inv = get_M_TGM_inv(A)

    u, iters = iterative_solve(A, b, zeros(N - 1), M_TGM_inv, 1e-6, 1000)

    # Plot solution
    u_exact(x) = exp(x) * (1 - x)
    x = range(1 / N, 1 - 1 / N, length=N - 1)

    p = plot(x, u, label="Numerical Solution", title="TGM iterative Solution", xlabel="x", ylabel="u(x)")
    plot!(x, u_exact.(x), label="Exact Solution")
end

# Task 2
plot_absolute_values_combined(get_B_TGM, "\n" * L"Eigenvalues of $B_{TGM}$" * "\n")
savefig("plots/task_12_eigenvalues.png")


# Task 3
plot_spectral_heatmap(get_B_TGM)
savefig("plots/task_12_spectral_heatmap.png")

# Task 4
# plot_convergence_combined(c_values, N_values, get_B_TGM, "Multigrid", get_M_TGM_inv)
plot_convergece_rate_against_spectral_radius(c_values, N_values, get_B_TGM, "Multigrid", get_M_TGM_inv)
savefig("plots/task_12_convergence.png")

