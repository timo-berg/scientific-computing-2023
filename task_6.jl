using LinearAlgebra
using Plots
using ColorSchemes

# TODO: Looks generally fine but still a bit fishy. Why does the solution not reach the right boundary?
# TODO: Tolerance is reached but error is still high. Why?

include("utils.jl")
include("gauss_seidel.jl")
include("conjugate_gradient.jl")
include("plot_def.jl")

function simulate(N, c, get_M_inv, tol=1e-10, max_iter=10000)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u0 = zeros(N - 1)
    M_inv = get_M_inv(A)
    # u, err, iter, res = conjugate_gradient(M_inv * A, M_inv * b, u0, tol, max_iter)
    u, err, iter, res = preconditioned_cg(A, b, u0, tol, max_iter, M_inv)
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
    x_numeric = range(1 / N, 1 - 1 / N, length=N - 1)
    p = plot(x_numeric, u, label="Numerical Solution", title="Conjugate Gradient Solution", xlabel="x", ylabel="u(x)")
    plot!(x_numeric, u_exact.(x_numeric), label="Exact Solution")
    print("Error: ", norm(u_exact.(x_numeric) - u), "\n")
    print("Iterations: ", iter)
    return p
end

function plot_convergences(get_M_inv, method)
    # Plot convergence for different values of N
    p_N = plot(xlabel="Iteration", ylabel="Log Error")
    N_values = [100, 200, 300, 1000]
    c = 60
    for N in N_values
        u, err, iter = simulate(N, c, get_M_inv)
        plot!(p_N, err, label="N = $N, c = $c", yscale=:log10)
    end

    # Plot convergence for different values of c
    p_c = plot(xlabel="Iteration", ylabel="Log Error")
    c_values = [0, 20, 40, 60, 80, 100, 120]
    N = 100
    for c in c_values
        u, err, iter = simulate(N, c, get_M_inv)
        # Use viridis colors for C
        color_palette = ColorSchemes.bamako
        color = color_palette.colors[Int(round((c - 0) / (120 - 0) * (length(color_palette.colors) - 1)) + 1)]
        plot!(p_c, err, label="N = $N, c = $c", yscale=:log10, legend=:bottomleft, color=color)
    end

    # Combined plot
    title_plot = plot(title="Convergence of Conjugate Gradient for $method", grid=false, showaxis=false, bottom_margin=-30Plots.px)
    p = plot(title_plot, p_N, p_c, layout=@layout([A{0.01h}; B; C]), size=(800, 800))

    return p
end
# plot_solution(100, 10)
plot_convergences(get_inv_M_SGS, "Gauss-Seidel")

savefig("plots/task_6_convergence.png")