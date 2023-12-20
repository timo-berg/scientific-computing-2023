using LinearAlgebra
using Plots
using LaTeXStrings

###### Works but not sure if its correct? 

N = 300
c = 60
h = N - 1
H = Int(N / 2 - 1)


B_CGC = get_B_CGC(H, h, c)

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
plot(p, q, layout=(1, 2), size=(1000, 500))