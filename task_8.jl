using LinearAlgebra
using Plots
using LaTeXStrings
include("multi_grid.jl")

###### Works but not sure if its correct? 

N = 300 # Nr of mesh elements. 
c = 60
A = get_A(N, c)
# n_h = N - 1 # Nr of internal nodes in fine grid
# n_H = Int(N / 2 - 1) # Nr of exernal nodes in coarse grid

B_CGC = get_B_CGC(A)

B_eigvals = eigvals(B_CGC)

# Absoulte values of eigenvalues in a plot
p = plot(title=L"Absolute Eigenvalues of $B_{CGC}$")
abs_eigvals = map(x -> abs(x), B_eigvals)
scatter!(p, abs_eigvals, label="", ms=2, markerstrokewidth=0)
p

