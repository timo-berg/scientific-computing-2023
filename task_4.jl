# Run Gauss-Seidel as a solver. Determine the asymptotic rate of convergence. Compare
# with the spectral radius.
# TODO: Get the spectral radius from task 3 and plot the convergence rate 
# for different values of c and N
# TODO: error gets worse at some

using LinearAlgebra
using Plots
# using Debugger
include("utils.jl")
include("task_3.jl")
include("plot_def.jl")

function do_step(u, A, b)
    u_old = copy(u)
    n = length(u)

    for i in 1:n
        sum1 = i > 1 ? sum(A[i, j] * u[j] for j in 1:i-1) : 0
        sum2 = i != n ? sum(A[i, j] * u_old[j] for j in i+1:n) : 0
        u[i] = (b[i] - sum1 - sum2) / A[i, i]
    end

    return u
end

function gauss_seidel(A, b, u0, tol=1e-6, max_iter=1000)
    u = u0
    err_list = []
    u_new = do_step(u, A, b)
    err = norm(A * u_new - b)
    push!(err_list, err)
    iter = 1
    while err > tol && iter < max_iter
        u = u_new
        u_new = do_step(u, A, b)

        err = norm(A * u_new - b) / norm(b)

        push!(err_list, err)
        iter += 1
    end

    return u_new, err_list, iter
end

function simulate_gs(N, c, tol=1e-10, max_iter=10000)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u0 = zeros(N - 1)
    u, err, iter = gauss_seidel(A, b, u0, tol, max_iter)
    return u, err, iter
end

function simulate_iterative(N, c, get_M_inv, tol=1e-10, max_iter=10000)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u0 = zeros(N - 1)
    M_inv = get_M_inv(A)
    # The +1 is required because the GS has max_iter + 1 iterations implicitly
    u, iter, residuals = iterative_solve(A, b, u0, M_inv, tol, max_iter + 1, true)

    err = map(x -> norm(x) / norm(b), residuals)

    return u, err, iter
end

function plot_convergence(N, c, tol=1e-6, max_iter=1000)
    u, err, iter = simulate(N, c, tol, max_iter)
    p = plot(err, title="Convergence of Gauss-Seidel", xlabel="Iteration", ylabel="Error", label="N = $N, c = $c", yscale=:log10)
    return p
end

# Plots the solution for multiple specified iterations 
function plot_multiple_sols(N,c, tol=1e-6)
    iterations = [10,50,100,135,200, 270]
    u_exact(x) = exp(x) * (1 - x)
    x_numeric = range(1 / N, 1 - 1 / N, length=N - 1)
    p = plot(label="Numerical Solutions", title="Gauss Seidel Solutions vs. Numerical Solution", xlabel="x", ylabel="u(x)")
    plot!(x_numeric, u_exact.(x_numeric), label="Exact Solution")
    # Color gradinet for the iterations
    colors = cgrad(:blues, length(iterations), rev=true)
    for (i, iteration) in enumerate(reverse(iterations))
        u, err, iter = simulate_gs(N, c, tol, iteration)
        # Plot the u values with color depnedning on the iteration number
        plot!(p, x_numeric, u, label="Iteration: $iteration", color=colors[i])
    end
    return p
end

function plot_solution(N, c, tol=1e-6, max_iter=80)
    u_exact(x) = exp(x) * (1 - x)
    u, err, iter = simulate_gs(N, c, tol, max_iter)
    x_numeric = range(1 / N, 1 - 1 / N, length=N - 1)
    p = plot(x_numeric, u, label="Numerical Solution", title="Exact Solution vs. Numerical Solution", xlabel="x", ylabel="u(x)")
    plot!(x_numeric, u_exact.(x_numeric), label="Exact Solution")
    print("Error: ", norm(u_exact.(x_numeric) - u), "\n")
    print("Iterations: ", iter)
    return p
end

function plot_convergence_combined(get_B, method, get_M_inv=nothing, tol=1e-6, max_iter=1000)
    c_values = 0:25:100
    N_values = [50]

    p = plot(title="Rate of convergence for $method", ylabel="Error", xlabel="Iteration")

    spectral_values = get_spectral_Bgs_values(c_values, N_values, get_B)

    for (i_N, N) in enumerate(N_values)
        for (i_c, c) in enumerate(c_values)
            if method == "Gauss-Seidel"
                u, err, iter = simulate_gs(N, c, tol, max_iter)
            else
                u, err, iter = simulate_iterative(N, c, get_M_inv, tol, max_iter)
            end
            spectral_value = spectral_values[i_c, i_N]
            plot!(p, err, label="N = $N, c = $c", yscale=:log10)

            # Annotate the spectral radius
            # Round the spectral value to 3 decimals
            spectral_value = round(spectral_value, digits=3)
            annotate!(p, iter, err[end], text("$spectral_value", 7, :left))
        end
    end
    return p
end


function plot_convergece_rate_against_spectral_radius(get_B, method, M_inv=I, tol=1e-6, max_iter=1000)
    c_values = 0:10:100
    N_values = [20, 50]

    spectral_values = get_spectral_Bgs_values(c_values, N_values, get_B)
    convergence_rates = zeros(length(c_values), length(N_values))

    p = plot(title="Log Convergence Rate vs. Spectral Radius", xlabel="Spectral Radius", ylabel="Log Convergence Rate")
    for (i_N, N) in enumerate(N_values)
        for (i_c, c) in enumerate(c_values)
            if method == "Gauss-Seidel"
                u, err, iter = simulate_gs(N, c, tol, max_iter)
            else
                u, err, iter = simulate_iterative(N, c, M_inv, tol, max_iter)
            end

            # Determine the log slope of the last 10 iterations
            log_err = log10.(err[end-10:end])
            convergence_rates[i_c, i_N] = (log_err[end] - log_err[1]) / 10

        end
    end

    # Plot the convergence rates against the spectral radius as a scatter plot
    for (i_N, N) in enumerate(N_values)
        scatter!(p, spectral_values[:, i_N], convergence_rates[:, i_N], label="N = $N", markerstrokewidth=0)
        # Annotate each point with the corresponding c value
        for (i_c, c) in enumerate(c_values)
            annotate!(p, spectral_values[i_c, i_N]+0.004, convergence_rates[i_c, i_N]-0.002, text("$c", 7, :left))
        end
    end

    return p
end

# plot_multiple_sols(100, 100)

# plot_convergence_combined(get_BGS, "Gauss-Seidel")


plot_convergece_rate_against_spectral_radius(get_BGS, "Gauss-Seidel")

# savefig("plots/task_4_error_vs_iterations_N100.png")
# savefig("plots/task_4_convergence_rate.png")
# savefig("plots/task_4_solutions_sample.png")