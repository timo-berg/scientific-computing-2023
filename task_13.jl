include("multi_grid.jl")
include("task_5.jl")
include("task_6.jl")
include("task_7.jl")
include("plot_def.jl")

# Task 5
compare_condition_numbers(get_M_TGM_inv, "TGM")
savefig("plots/task_13_condition_numbers.png")

# Task 6
plot_convergences(get_M_TGM_inv, "Multigrid")
savefig("plots/task_13_convergence.png")

# Task 7
plot_ritz_values(get_M_TGM_inv, "Multigrid")
