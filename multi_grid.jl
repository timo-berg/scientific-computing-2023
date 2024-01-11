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

function get_fine_and_coarse_nr_node(N)
    n_h, n_H = N - 1, Int(N / 2 - 1)
    return n_h, n_H
end

function get_M_CGC_inv(A)
    n_h, n_H = get_fine_and_coarse_nr_node(size(A, 1) + 1)
    I_H_to_h = linear_interpolation(n_h, n_H)
    I_h_to_H = half_weight_mapping(n_h, n_H)
    A_H = I_h_to_H * A * I_H_to_h

    return I_H_to_h * inv(A_H) * I_h_to_H
end

function get_B_CGC(A)
    M_CGC_inv = get_M_CGC_inv(A)

    return I - M_CGC_inv * A
end


function get_B_TGM(A)
    B_CGC = get_B_CGC(A)
    B_BGS = get_B_BGS(A)
    B_FGS = get_B_FGS(A)

    return B_BGS * B_CGC * B_FGS
end

function get_M_TGM_inv(A)
    M_BGS = get_inv_M_BGS(A)
    M_CGC = get_M_CGC_inv(A)
    M_FGS = get_inv_M_FGS(A)

    return M_BGS + M_CGC + M_FGS - M_BGS * A * M_CGC - M_CGC * A * M_FGS - M_BGS * A * M_FGS + M_BGS * A * M_CGC * A * M_FGS
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


function test_TGM(A)
    B_TGM = get_B_TGM(A)
    M_TGM_inv = get_M_TGM_inv(A)

    return norm(B_TGM - I + M_TGM_inv * A)
end
