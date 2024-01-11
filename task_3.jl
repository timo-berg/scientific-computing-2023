include("gauss_seidel.jl")
include("plot_def.jl")

function get_spectral_Bgs_values(c_values, n_values, get_B)
    num_c = length(c_values)
    num_n = length(n_values)
    spectral_radiuses = Matrix{Float64}(undef, num_c, num_n)

    for c_i = 1:num_c
        for n_i = 1:num_n
            A = get_A(n_values[n_i], c_values[c_i])
            eigenvalues = eigvals(get_B(A))
            abs_eig = map(x -> abs(x), eigenvalues)
            spectral_radius = maximum(abs_eig)
            spectral_radiuses[c_i, n_i] = spectral_radius
        end
    end
    return spectral_radiuses
end

function plot_spectral_heatmap(get_B)
    c_values = [-5, 0, 5, 10, 100]
    n_values = [10, 50, 100, 200, 300]
    spectral_radiuses = get_spectral_Bgs_values(c_values, n_values, get_B)

    c_indices = 1:length(c_values)
    n_indices = 1:length(n_values)
    p = heatmap(n_indices, c_indices, spectral_radiuses, xticks=(n_indices, n_values), yticks=(c_indices, c_values))

    title!(p, "Spectral Radius Values")
    xlabel!(p, "N")
    ylabel!(p, "c")

    ann = [(j, i, text(round(spectral_radiuses[i, j], digits=4), 10, :white, :center))
           for i in 1:length(c_values) for j in 1:length(n_values)]
    annotate!(ann, linecolor=:white)
end

plot_spectral_heatmap(get_BGS)

savefig("plots/task_3_spectral_radius_heatmap.png")