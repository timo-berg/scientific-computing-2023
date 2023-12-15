# Run Gauss-Seidel as a solver. Determine the asymptotic rate of convergence. Compare
# with the spectral radius.
# TODO: Get the spectral radius from task 3 and plot the convergence rate 
# for different values of c and N

using LinearAlgebra
using Plots
include("utils.jl")

function do_step(u, A, b)
    u_old = copy(u)
    n = length(u)

    for i in 1:n
        sum1 = i > 1 ? sum(A[i, j] * u[j] for j in 1:i-1) : 0
        sum2 = i != n ? sum(A[i, j] * u_old[j] for j in i+1:n) : 0
        u[i] = (b[i] - sum1 - sum2) / A[i, i]
    end

    err = maximum(abs.(u_old - u))

    return u, err
end

function gauss_seidel(A, b, u0, tol=1e-6, max_iter=1000)
    u = u0
    err_list = []
    u_new, err = do_step(u, A, b)
    iter = 1
    while err > tol && iter < max_iter
        u = u_new
        u_new, err = do_step(u, A, b)
        push!(err_list, err)
        iter += 1
    end
    return u_new, err_list, iter
end

function simulate(N, c, tol=1e-10, max_iter=10000)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u0 = zeros(N - 1)
    u, err, iter = gauss_seidel(A, b, u0, tol, max_iter)
    return u, err, iter
end

function plot_convergence(N, c, tol=1e-6, max_iter=1000)
    u, err, iter = simulate(N, c, tol, max_iter)
    p = plot(err, title="Convergence of Gauss-Seidel", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c")
    return p
end

function plot_solution(N, c, tol=1e-10, max_iter=100000)
    u_exact(x) = exp(x) * (1 - x)
    u, err, iter = simulate(N, c, tol, max_iter)
    x_numeric = range(0, 1, length=N - 1)
    p = plot(x_numeric, u, label="Numerical Solution", title="Exact Solution vs. Numerical Solution", xlabel="x", ylabel="u(x)")
    plot!(x_numeric, u_exact.(x_numeric), label="Exact Solution")
    print("Error: ", norm(u_exact.(x_numeric) - u), "\n")
    print("Iterations: ", iter)
    return p
end

# plot_convergence(100, 2)
plot_solution(100, 2)

