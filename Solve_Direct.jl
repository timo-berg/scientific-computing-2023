using LinearSolve
using Plots
include("utils.jl")


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