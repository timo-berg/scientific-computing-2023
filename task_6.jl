using LinearAlgebra
using Plots
using ColorSchemes

include("utils.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")

function simulate(N, c, tol=1e-10, max_iter=10000)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u0 = zeros(N - 1)
    M_inv = get_inv_M_SGS(A)
    u, err, iter = conjugate_gradient(A, b, u0, tol, max_iter)
    # u, err, iter = preconditioned_cg(A, b, u0, tol, max_iter, M_inv)
    return u, err, iter
end

function plot_convergence(N, c, tol=1e-6, max_iter=1000)
    u, err, iter = simulate(N, c, tol, max_iter)
    p = plot(err, title="Convergence of Conjugate Gradient", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c")
    return p
end

function plot_solution(N, c, tol=1e-10, max_iter=100000)
    u_exact(x) = exp(x) * (1 - x)
    u, err, iter = simulate(N, c, tol, max_iter)
    x_numeric = range(0, 1, length=N - 1)
    p = plot(x_numeric, u, label="Numerical Solution", title="Conjugate Gradient Solution", xlabel="x", ylabel="u(x)")
    plot!(x_numeric, u_exact.(x_numeric), label="Exact Solution")
    print("Error: ", norm(u_exact.(x_numeric) - u), "\n")
    print("Iterations: ", iter)
    return p
end


# Plot convergence for different values of N
p_N = plot()
N_values = [100, 200, 300, 1000]
c = 60
for N in N_values
    u, err, iter = simulate(N, c)
    plot!(p_N, err, label="N = $N, c = $c", yscale=:log10)
end

# Plot convergence for different values of c
p_c = plot()
c_values = [0, 20, 40, 60, 80, 100, 120]
N = 100
for c in c_values
    u, err, iter = simulate(N, c)
    # Use viridis colors for C
    color_palette = ColorSchemes.broc
    color = color_palette.colors[Int(round((c - 0) / (120 - 0) * (length(color_palette.colors) - 1)) + 1)]
    plot!(p_c, err, label="N = $N, c = $c", yscale=:log10, legend=:bottomleft, color=color)
end

# Combined plot
plot(p_N, p_c, layout=(2, 1), size=(800, 800))