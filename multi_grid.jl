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

function get_B_CGC(H, h, c)
    I_H_to_h = linear_interpolation(h, H)
    I_h_to_H = half_weight_mapping(h, H)
    A = get_A(h + 1, c)

    A_H_inv = inv(I_h_to_H * A * I_H_to_h)

    return I - I_H_to_h * A_H_inv * I_h_to_H * A
end