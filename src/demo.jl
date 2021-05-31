proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"Utility.jl")
include(proj_dir*"cpMatrices.jl")
include(proj_dir*"Moments.jl")
include(proj_dir*"Constraints.jl")
include(proj_dir*"cpModel.jl")
include(proj_dir*"Compute.jl")
# include(proj_dir*"Batch_proc.jl")

using .Utility
using .cpMatrices
using .Moments
using .Constraints
using .cpModel
using .Compute


cp_mats = cpMatrices.get_cp_mats() # load the matrices
A = cp_mats["M6"]                  # Pick a specific one
t = 2                              # Choose the level of the hierarchy
conlist = "DagXXwGsG"              # Choose which additional constraints to add
mod = cpModel.Modelξₜᶜᵖ(A,t,conlist) # Build the model
Compute.Computeξₜᶜᵖ(mod)  # Solve the model using Mosek
