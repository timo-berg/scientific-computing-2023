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

function plot_histogram_c_range()
    N = 100
    p = plot()
    c_values = [-5, -2, 0, 2, 5] # [-5, -1, 0, 1, 5]
    for c = c_values
        A = get_A(N, c)
        scatter!(p, eigvals(get_BGS(A)), label="C = $c", ms=3)
    end
    p
end

function get_inv_M_SGS(A)
    D = Diagonal(A)
    E = -LowerTriangular(A)
    F = -UpperTriangular(A)
    return inv(D - F) * D * inv(D - E)
end


function plot_histogram_n_range()
    c_value = 2
    p = plot()
    title!(p, L"Eigenvalues of $B_{GS}$ with $c = %$c_value$")
    n_values = [10, 100, 200, 400] # [-5, -1, 0, 1, 5]
    for n = n_values
        A = get_A(n, c_value)
        scatter!(p, eigvals(get_BGS(A)), label="N = $n", ms=3)
    end
    p
end

function plot_aboslute_values_c_range()
    n_value = 100
    c_values = [-5, 0, 5, 10, 100]
    p = plot()
    for c = c_values
        A = get_A(n_value, c)
        eigenvalues = eigvals(get_BGS(A))
        abs_eig = map(x -> abs(x), eigenvalues)
        scatter!(p, sort(abs_eig, rev=true), label=L"c=%$c", ms=3)
    end
    savefig(p, "plots/plot.pdf")
    p
end

function plot_aboslute_values_n_range()
    c_value = 10
    n_values = [10, 100, 200, 400]
    p = plot()
    for n = n_values
        A = get_A(n, c_value)
        eigenvalues = eigvals(get_BGS(A))
        abs_eig = map(x -> abs(x), eigenvalues)
        scatter!(p, sort(abs_eig, rev=true), label=L"N=%$n", ms=3)
    end
    savefig(p, "plots/plot.pdf")
    p
end

function test_get_BGS()
    A = [2 -1 0; -1 2 -1; 0 -1 2]
    B_expected = [0 1/2 0; 0 1/4 1/2; 0 1/8 1/4]
    print(get_BGS(A), B_expected)
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

plot_aboslute_values_n_range()


# plot_histogram_n_range()



