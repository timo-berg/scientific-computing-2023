using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")

###### Works but not sure if its correct? 

N = 300 # Nr of mesh elements. 
c = 60
n_h = N - 1 # Nr of internal nodes in fine grid
n_H = Int(N / 2 - 1) # Nr of exernal nodes in coarse grid

B_CGC = get_B_CGC(n_H, n_h , c)

B_eigvals = eigvals(B_CGC)

# Plot the eigenvalues with real part < 0.5
p = plot(title=L"Small Eigenvalues of $B_{CGC}$")
small_eigvals = [eigval for eigval in B_eigvals if real(eigval) < 0.5]
scatter!(p, real.(small_eigvals), imag.(small_eigvals), label="", ms=3)

# Plot the eigenvalues with real part > 0.5
q = plot(title=L"Large Eigenvalues of $B_{CGC}$")
large_eigvals = [eigval for eigval in B_eigvals if real(eigval) > 0.5]
scatter!(q, real.(large_eigvals), imag.(large_eigvals), label="", ms=3)


# Combine the plots
# plot(p, q, layout=(1, 2), size=(1000, 500))


# Absoulte values of eigenvalues in a plot
p = plot(title=L"Absolute Eigenvalues of $B_{CGC}$")
abs_eigvals = map(x -> abs(x), B_eigvals)
scatter!(p, sort(abs_eigvals, rev=true), label="", ms=3)
p