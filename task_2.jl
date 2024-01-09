using LinearAlgebra
using Plots


function plot_eig_complex_c_range(get_B)
    N = 100
    p = plot()
    c_values =  [0, 5, 10, 100]
    for c = c_values
        A = get_A(N, c)
        scatter!(p, eigvals(get_B(A)), label="C = $c", ms=3)
    end
    p
end


function plot_eig_complex_n_range(get_B)
    c_value = 0
    p = plot()
    title!(p, L"Eigenvalues of $B_{GS}$ with $c = %$c_value$")
    n_values = [10, 100, 200, 400] # [-5, -1, 0, 1, 5]
    for n = n_values
        A = get_A(n, c_value)
        scatter!(p, eigvals(get_B(A)), label="N = $n", ms=3)
    end
    p
end


function plot_aboslute_values_c_range(get_B)
    n_value = 100
    c_values = [0, 5, 10, 100]
    p = plot()
    for c = c_values
        A = get_A(n_value, c)
        eigenvalues = eigvals(get_B(A))
        abs_eig = abs.(eigenvalues)
        scatter!(p, sort(abs_eig, rev=true), label=L"c=%$c", ms=2, markerstrokewidth=0)
    end
    savefig(p, "plots/plot.pdf")
    p
end

function plot_aboslute_values_n_range(get_B)
    c_value = 0
    n_values = [10, 100, 200, 400]
    p = plot()
    for n = n_values
        A = get_A(n, c_value)
        eigenvalues = eigvals(get_B(A))
        abs_eig = map(x -> abs(x), eigenvalues)
        scatter!(p, sort(abs_eig, rev=true), label=L"N=%$n", ms=2, markerstrokewidth=0)
    end
    savefig(p, "plots/plot.pdf")
    p
end

function plot_absolute_values_combined(get_B, title)
    p1 = plot_aboslute_values_c_range(get_B)
    p2 = plot_aboslute_values_n_range(get_B)

    # Combine plots and set global title
    title = plot(title = title, grid = false, showaxis = false, bottom_margin = -30Plots.px)

    p = plot(title, p1, p2, layout = @layout([A{0.01h}; B; C]))
end

plot_absolute_values_combined(get_BGS, L"Eigenvalues of $B_{GS}$")

