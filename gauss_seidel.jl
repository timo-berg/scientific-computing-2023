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
println(get_A(3, 0))
println(get_BGS(get_A(3,0)))
eigvals(get_BGS(get_A(10,0)))

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
    A_lower = LowerTriangular(A)
    A_upper = UpperTriangular(A)
    return inv(A_upper) * D * inv(A_lower)
end


function plot_histogram_n_range()
    c_value = 0
    p = plot()
    title!(p, L"Eigenvalues of $B_{GS}$ with $c = %$c_value$")
    n_values = [10] #, 100, 200, 400] # [-5, -1, 0, 1, 5]
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
    c_value = 0
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

function get_spectral_Bgs_values(c_values, n_values)
    num_c = length(c_values)
    num_n = length(n_values)
    spectral_radiuses = Matrix{Float64}(undef, num_c, num_n)
    
    for c_i = 1:num_c
        for n_i = 1:num_n
        A = get_A(n_values[n_i], c_values[c_i])
        eigenvalues = eigvals(get_BGS(A))
        abs_eig = map(x -> abs(x), eigenvalues)
        spectral_radius = maximum(abs_eig)
        spectral_radiuses[c_i, n_i] = spectral_radius
        end
    end
    return spectral_radiuses
end

function plot_spectral_heatmap()
    c_values = -100:10:20
    n_values = StepRange{Int}(10, 5, 60)
    spectral_radiuses = get_spectral_Bgs_values(c_values, n_values)
    p = heatmap(n_values, c_values, spectral_radiuses)
    title!(p, "Spectral Radius Values")
    xlabel!(p, "N values")
    ylabel!(p, "C values")
end

plot_spectral_heatmap()

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

plot_aboslute_values_n_range()

#plot_histogram_n_range()