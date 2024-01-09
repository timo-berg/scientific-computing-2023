using LinearAlgebra
using Plots


function plot_eig_complex_c_range()
    N = 100
    p = plot()
    c_values = [-5, -2, 0, 2, 5] # [-5, -1, 0, 1, 5]
    for c = c_values
        A = get_A(N, c)
        scatter!(p, eigvals(get_BGS(A)), label="C = $c", ms=3)
    end
    p
end


function plot_eig_complex_n_range()
    c_value = 0
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
        abs_eig = abs.(eigenvalues)
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

plot_eig_complex_n_range()
# plot_aboslute_values_c_range()
# plot_aboslute_values_n_range()
