"""
    conjugate_gradient(A, b, x0, tol, maxiter)

Solves the linear system Ax = b using the Conjugate Gradient method.

Arguments:
- `A`: The coefficient matrix of the linear system.
- `b`: The right-hand side vector.
- `x0`: The initial guess for the solution.
- `tol`: The tolerance for convergence.
- `maxiter`: The maximum number of iterations.

Returns:
- `x`: The solution vector.
- `residuals`: The residuals at each iteration.
- `iter`: The number of iterations used.
"""
function conjugate_gradient(A, b, x0, tol, maxiter)
    # Standard implementation of the conjugate gradient method
    r = b
    r_old = p = 0 # Initialize to zero, will only be used in iteration >2
    if norm(r) / norm(b) < tol
        return x0
    end

    x = x0
    convergence = [norm(r)]
    residuals = [r]
    b_norm = norm(b)

    for i = 1:maxiter
        if i == 1
            p = r
        else
            β = dot(r, r) / dot(r_old, r_old)
            p = r + β * p
        end

        A_p = A * p
        α = dot(r, r) / dot(p, A_p)
        x = x + α * p
        r_old = r
        r = r - α * A_p

        push!(residuals, r)
        push!(convergence, norm(r) / b_norm)
        if norm(r) / norm(b) < tol
            real_residual = A * x - b
            if (norm(real_residual) / b_norm) > tol
                r = real_residual
            else
                return x, convergence, i, residuals
            end
        end
    end

    return x, convergence, maxiter, residuals
end

""" 
Preconditioned implementation of the conjugate gradient method
"""
function preconditioned_cg(A, b, x0, tol, maxiter, M_inv)
    r = b
    p = z_old = r_old = 0 # Initialize to zero, will only be used in iteration >2

    if norm(r) / norm(b) < tol
        return x0
    end

    x = x0
    convergence = [norm(r)]
    residuals = [r]

    for i = 1:maxiter
        z = M_inv * r
        if i == 1
            p = z
        else
            β = dot(r, z) / dot(r_old, z_old)
            p = z + β * p
        end

        # Updating the solution
        A_p = A * p
        α = dot(r, z) / dot(p, A_p)
        x = x + α * p
        z_old = z

        # Updating the residual
        r_old = r
        r = r - α * A_p
        push!(convergence, norm(r) / norm(b))
        push!(residuals, r)
        if norm(r) / norm(b) < tol
            return x, convergence, i, residuals
        end
    end

    return x, convergence, maxiter, residuals
end

function conjugate_gradient_solver(A, b)
    x0 = zeros(length(b))
    conjugate_gradient(A, b, x0, 1e-6, 1000)
end

function preconditioned_cg_solver(A, b, M_inv)
    x0 = zeros(length(b))
    preconditioned_cg(A, b, x0, 1e-6, 1000, M_inv)
end

#@time conjugate_gradient_solver(M_inv * A, M_inv * b)
#@time preconditioned_cg_solver(A, b, M_inv)