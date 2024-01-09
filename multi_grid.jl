include("utils.jl")
using LinearAlgebra

function linear_interpolation(h, H)
    I = zeros(h, H)

    for i = 2:H-1
        I[2*i-2, i] = 1
        I[2*i-1, i] = 2
        I[2*i, i] = 1
    end
    # When interpolating bedween the 'edge' inner points and the boundary value, we just choose the inner point.
    I[1, 1] = 2 
    I[2, 1] = 1
    I[h-1, H] = 1
    I[h, H] = 2

    return I * 1 / 2
end

function half_weight_mapping(h, H)
    I = zeros(H, h)
    # Interior points
    for i = 2:H-1
        I[i, 2*i-2] = 1
        I[i, 2*i-1] = 2
        I[i, 2*i] = 1
    end

    # Boundary points
    I[1, 1] = 4
    I[H, h] = 4

    return I * 1 / 4
end

print(half_weight_mapping(5, 2))


"""
Returns the number of nodes in the 1 dimenionsal fine and double coarse grid given the number of meshes.
"""
function get_fine_and_coarse_nr_node(N)
    n_h, n_H = N-1, Int(N/2)
    return n_h, n_H
end

"""
Returns the inverse of the coarse grid correction operator
"""
function get_M_CGC_inv(n_H, n_h, c)
    I_H_to_h = linear_interpolation(n_h, n_H)
    I_h_to_H = half_weight_mapping(n_h, n_H)
    A_h = get_A(n_h + 1, c) # +1 because function takes nr of meshes, not nr of internal points
    A_H = I_h_to_H * A_h * I_H_to_h

    return I_H_to_h * inv(A_H) * I_h_to_H
end

"""
    get_B_CGC(n_H, n_h, c)

Returns the error propagation matrix for the CGC method.

# Arguments
- `n_H::Int`: Number of internal nodes in the coarse grid
- `n_h::Int`: Number of internal nodes in the fine grid
- `c::Number`: The c value in the diff equation
- `A::Matrix`: The coefficient matrix for the diff equation
"""
function get_B_CGC(n_H, n_h, c, A)
    M_CGC_inv = get_M_CGC_inv(n_H, n_h, c)

    return I - M_CGC_inv * A
end

# Two-grid cycle
function two_grid_cycle(A, b, u, I_H_to_h, I_h_to_H, preconditioner, solver=direct_solve)
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
    e_H = solver(A_H, r_H)

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
        u = two_grid_cycle(A, b, u, I_H_to_h, I_h_to_H, M_inv, solver)

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
