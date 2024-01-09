include("gauss_seidel.jl")

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
    c_values = 0:10:100
    n_values = StepRange{Int}(10, 10, 60)
    spectral_radiuses = get_spectral_Bgs_values(c_values, n_values, get_B)
    p = heatmap(n_values, c_values, spectral_radiuses)
    title!(p, "Spectral Radius Values")
    xlabel!(p, "N values")
    ylabel!(p, "C values")
end

plot_spectral_heatmap(get_BGS)