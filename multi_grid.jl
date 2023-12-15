include("utils.jl")


function linear_interpolation(H, h)
    I = zeros(H, h)
    # First column 
    I[1, 1] = 2
    I[2, 1] = 1


    # Interior points
    for i = 2:h
        I[2*i-2, i] = 1
        I[2*i-1, i] = 2
        I[2*i, i] = 1
    end

    return I * 1 / 2
end

function half_weight_mapping(H, h)
    I = zeros(h, H)
    # First row 
    I[1, 1] = 2
    I[1, 2] = 1


    # Interior points
    for i = 2:h
        I[i, 2*i-2] = 1
        I[i, 2*i-1] = 2
        I[i, 2*i] = 1
    end

    return I * 1 / 4
end

function get_B_CGC(H, h, c)
    I_h_to_H = linear_interpolation(H, h)
    I_H_to_h = half_weight_mapping(H, h)
    A = get_A(h + 1, c)

    A_H_inv = inv(I_h_to_H * A * I_H_to_h)

    return I - I_H_to_h * A_H_inv * I_h_to_H * A
end
