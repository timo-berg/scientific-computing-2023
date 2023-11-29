using Plots
include("utils.jl")


# Set the coefficients
α = 1
β = 0
c = -1
f(x) = -2 * exp(x)


# Correct solution
u_exact(x) = exp(x) * (1 - x)

mesh_sizes = [3, 6, 10]

### Plot exact solution
x_exact = range(0, 1, length=1000)
solution_plot = plot(x_exact, u_exact.(x_exact), label="Exact Solution", title="Exact Solution vs. Numerical Solution", xlabel="x", ylabel="u(x)")

# Vary the mesh size
for N in mesh_sizes
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


mesh_sizes = 3:50

### Investigate the error
errors = []
for N in mesh_sizes
    # Get the coefficient matrix A and the right hand side vector b
    A = get_A(N, c)
    b = get_b(N, f, α, β)

    # Solve the linear system Ax = b
    u = A \ b

    # Plot the solution
    x_numeric = range(0, 1, length=N + 1)

    # Add error to the list
    push!(errors, maximum(abs.(u - u_exact.(x_numeric[2:end-1]))))
end

# Plot the errors
error_plot = scatter(mesh_sizes, errors, title="Discretization Error scaling with number of mesh points", xlabel="Mesh Points", ylabel="Error", legend=false)
# plot!(mesh_sizes, 1 ./ mesh_sizes.^4)

display(solution_plot)
display(error_plot)