include("utils.jl")
using LinearAlgebra

function linear_interpolation(h, H)
    I = zeros(h, H)

    for i = 1:H
        I[2*i-1, i] = 1
        I[2*i, i] = 2
        I[2*i+1, i] = 1
    end

    return I * 1 / 2
end

function half_weight_mapping(h, H)
    I = zeros(H, h)
    # Interior points
    for i = 1:H
        I[i, 2*i-1] = 1
        I[i, 2*i] = 2
        I[i, 2*i+1] = 1
    end

    return I * 1 / 4
end

"""
    get_B_CGC(n_H, n_h, c)

Returns the error propagation matrix for the CGC method.

# Arguments
- `n_H::Int`: Number of internal nodes in the coarse grid
- `n_h::Int`: Number of internal nodes in the fine grid
"""
function get_B_CGC(n_H, n_h, c)
    I_H_to_h = linear_interpolation(n_h, n_H)
    I_h_to_H = half_weight_mapping(n_h, n_H)
    A = get_A(n_h + 1, c)

    A_H_inv = inv(I_h_to_H * A * I_H_to_h)

    return I - I_H_to_h * A_H_inv * I_h_to_H * A
end

# Two-grid cycle
function two_grid_cycle(A, b, u, I_H_to_h, I_h_to_H, preconditioner)
    # Smoother
    S_h = I - preconditioner * A
    # Interpolated problem
    A_H = I_h_to_H * A * I_H_to_h
    
    # Pre smoothing
    u = S_h * u + preconditioner * b

    # residual
    r_h = b - A * u

    # Restrict residual
    r_H = I_h_to_H * r_h

    # Solve error equation
    e_H = zeros(length(r_H))
    e_H = inv(A_H) * r_H

    # Interpolate error
    e_h = I_H_to_h * e_H

    # Correct
    u = u + e_h

    # Post smoothing
    u = S_h * u + preconditioner * b

    return u
end


    
# Multi-grid cycle
function multigrid(A, b, ϵ, max_iter)
    n_h = length(b)
    n_H = Int((h + 1) / 2 - 1)

    # Get the projector
    I_H_to_h = linear_interpolation(n_h, n_H)
    I_h_to_H = half_weight_mapping(n_h, n_H)

    # Get the smoother
    M_inv = get_inv_M_SGS(A)

    # Initial guess
    u = zeros(h)

    errors = []
    
    # Iterate
    for i = 1:max_iter
        # Two-grid cycle
        u = two_grid_cycle(A, b, u, I_H_to_h, I_h_to_H, M_inv)

        # Error
        error = norm(A * u - b)
        push!(errors, error)
        # Check convergence
        if error < ϵ
            println("Converged after $i iterations")
            break
        end
    end

    println("Error: ", norm(A * u - b))

    return u, errors
end