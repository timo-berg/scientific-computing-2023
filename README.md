# scientific-computing-2023

When cloning the repository make sure you use
```
julia> using Pkg
julia> Pkg.activate("path/to/cloned_repository")
julia> Pkg.instantiate()
```
to use this repositories virtual environment.

If you add packages to the environment, use 
```
julia> Pkg.update()
```
and commit the `Project.toml` and `Manifest.toml` files.