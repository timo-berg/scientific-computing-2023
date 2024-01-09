include("multi_grid.jl")
include("task_5.jl")
include("task_6.jl")
include("task_7.jl")

# Task 5
compare_condition_numbers(get_M_TGM_inv, "TGM")

# Task 6
plot_convergences(get_M_TGM_inv, "Multigrid")

# Task 7
plot_ritz_values(get_M_TGM_inv, "Multigrid")
