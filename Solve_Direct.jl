using LinearSolve
using Plots

# Get the coefficient matrix A
function get_A(N, c)
    h = 1 / N
    A = zeros(N - 1, N - 1)

    # Set the boundary points
    A[1, 1] = c - 2 / h^2
    A[1, 2] = 1 / h^2

    A[N-1, N-2] = 1 / h^2
    A[N-1, N-1] = c - 2 / h^2



    # Set the interior points
    for i = 2:N-2
        A[i, i-1] = 1 / h^2
        A[i, i] = c - 2 / h^2
        A[i, i+1] = 1 / h^2
    end

    return A
end

# Get the right hand side vector b
function get_b(N, f, α=1, β=0)
    h = 1 / N

    b = zeros(N - 1)

    # Set the boundary points
    b[1] = f(h) - α / h^2
    b[N-1] = f(1 - h) - β / h^2

    # Set the interior points
    for i = 2:N-2
        b[i] = f(i * h)
    end

    return b
end

# Set the coefficients
α = 1
β = 0
c = -1
f(x) = -2 * exp(x)


# Correct solution
u_exact(x) = exp(x) * (1 - x)

### Plot exact solution
x_exact = range(0, 1, length=1000)
plot(x_exact, u_exact.(x_exact), label="Exact Solution")



# Vary the mesh size
for N in [3, 6, 10]
    # Get the coefficient matrix A and the right hand side vector b
    A = get_A(N, c)
    b = get_b(N, f, α, β)

    # Solve the linear system Ax = b
    u = A \ b

    # Plot the solution
    x_numeric = range(0, 1, length=N + 1)

    plot!(x_numeric, [α; u; β], label="N = $N")
end

plot!(legend=:topright)


### Investigate the error
for N in [3, 6, 10]
    # Get the coefficient matrix A and the right hand side vector b
    A = get_A(N, c)
    b = get_b(N, f, α, β)

    # Solve the linear system Ax = b
    u = A \ b

    # Plot the solution
    x_numeric = range(0, 1, length=N + 1)

    errors = u - u_exact.(x_numeric[2:end-1])

    println("N = $N, error = $(maximum(abs.(errors)))")
end