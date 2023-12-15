# Get the coefficient matrix A for the Poisson equation
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

# Returns the f function given a certain c
function construct_F(c)
    return function f(x)
        exp(x) * (c - 1) - x * exp(x) * (1 + c)
    end
end