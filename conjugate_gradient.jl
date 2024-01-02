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
    r = b - A * x0
    if norm(r) / norm(b) < tol
        return x0
    end

    p = r
    x = x0
    convergence = [norm(r)]

    for i = 1:maxiter
        A_p = A * p
        α = dot(r, r) / dot(p, A_p)
        x = x + α * p
        r_old = r
        r = r - α * A_p
        push!(convergence, norm(r))
        if norm(r) / norm(b) < tol
            return x, convergence, i
        end
        β = dot(r, r) / dot(r_old, r_old)
        p = r + β * p
    end

    return x, convergence, maxiter
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
        push!(convergence, norm(r))
        if norm(r) / norm(b) < tol
            return x, convergence, i
        end
    end

    return x, convergence, maxiter
end

function preconditioned_cg_residuals(A, b, x0, tol, maxiter, M_inv)
    # Preconditioned implementation of the conjugate gradient method
    r = b # With initial guess u=0 then the residual will be b.
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
        push!(convergence, norm(r))
        push!(residuals, r)
        if norm(r) / norm(b) < tol
            return residuals
        end
    end

    return residuals
end