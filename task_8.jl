using LinearAlgebra
using Plots

# TODO: gives error because inverse of A_H is singular. Probably error in multi_grid.jl
B_CGC = get_B_CGC(H, h, c)
