function iterative_solve(A, b, x_0, M_inv, tol=1e-6, maxiter=1000, return_residuals=false)
    r = A * x_0 - b
    x = x_0

    if return_residuals
        residuals = [r]
    end

    for i = 1:maxiter
        x = x - M_inv * r
        r = A * x - b
        if return_residuals
            push!(residuals, r)
        end
        if norm(r) <= tol
            if return_residuals
                return x, i, residuals
            else
                return x, i
            end
        end
    end

    if return_residuals
        return x, maxiter, residuals
    end
    return x, maxiter
end

# Get the coefficient matrix A for the Poisson equation
function get_A(N::Int, c::Number)
    h = 1 / N
    A = zeros(N - 1, N - 1)

    # Set the boundary points
    A[1, 1] = c + 2 / h^2
    A[1, 2] = -1 / h^2

    A[N-1, N-2] = -1 / h^2
    A[N-1, N-1] = c + 2 / h^2



    # Set the interior points
    for i = 2:N-2
        A[i, i-1] = -1 / h^2
        A[i, i] = c + 2 / h^2
        A[i, i+1] = -1 / h^2
    end

    return Symmetric(A)
end

# Get the right hand side vector b
function get_b(N, f, α=1, β=0)
    h = 1 / N

    b = zeros(N - 1)

    # Set the boundary points
    b[1] = f(h) + α / h^2
    b[N-1] = f(1 - h) + β / h^2

    # Set the interior points
    for i = 2:N-2
        b[i] = f(i * h)
    end

    return b
end

# Returns the f function given a certain c
function construct_F(c)
    return function f(x)
        exp(x) * (c + 1) + x * exp(x) * (1 - c)
    end
end


function check_correctness(N, c)
    A = get_A(N, c)
    b = get_b(N, construct_F(c))
    u_exact(x) = exp(x) * (1 - x)
    u = u_exact.(range(1 / N, 1 - 1 / N, length=N - 1))
    error = norm(A * u - b)

    return error
end

# """ Inspects the incfluence of c on the shape of F.
# """
# function inspect_c_influence()
#     N = 100
#     c_values = [0, 1, 10, 100]
#     p = plot(title="Influence of c on F", xlabel="x", ylabel="F(x)")
#     for c in c_values
#         F = construct_F(c)
#         x = range(0, 1, length=100)
#         plot!(p, x, F.(x), label="c = $c")
#     end
#     p
# end

# inspect_c_influence()