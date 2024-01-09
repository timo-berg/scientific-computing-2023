using LinearAlgebra
using Plots
using LaTeXStrings
include("utils.jl")

# Gets the Gauss Seidel error propagation matrix for a given A
function get_BGS(A)
    M = LowerTriangular(A)
    B_GS = I - inv(M) * A
    return B_GS
end

function get_B_BGS(A)
    M = UpperTriangular(A)
    B_GS = I - inv(M) * A
    return B_GS
end

function get_B_FGS(A)
    return get_BGS(A)
end

function get_inv_M_SGS(A)
    D = Diagonal(A)
    A_lower = LowerTriangular(A)
    A_upper = UpperTriangular(A)
    return inv(A_upper) * D * inv(A_lower)
end

function get_inv_M_FGS(A)
    M = LowerTriangular(A)
    return inv(M)
end

function get_inv_M_BGS(A)
    M = UpperTriangular(A)
    return inv(M)
end


function test_get_BGS()
    A = [2 -1 0; -1 2 -1; 0 -1 2]
    B_expected = [0 1/2 0; 0 1/4 1/2; 0 1/8 1/4]
    println(eigvals(get_BGS(A)))
    println(get_BGS(A),  B_expected)
end



# Test f at some points
function test_construct_F()
    f = construct_F(-1)
    expected_f(x) = -2 * exp(x)

    print(f(42) == expected_f(42))

    f_2 = construct_F(3)
    expected_f_2(x) = 2 * exp(x) - 4 * x * exp(x)

    print(f_2(42) == expected_f_2(42))

end
