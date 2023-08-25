This package is based on the work in the following publication:
> M. Korda, M. Laurent, V. Magron, and A. Steenkamp. 
> Exploiting ideal-sparsity in the generalized moment problem with application to matrix factorization ranks.
> Mathematical Programming, 2023. https://doi.org/10.1007/s10107-023-01993-x


This code models the problem of bounding the completely positive rank of a symmetric nonnegative matrix via moment hierarchies.

### Package requirements:
- [Test.jl](https://docs.julialang.org/en/v1/stdlib/Test/)
- [JuMP.jl](https://jump.dev/JuMP.jl/stable/)
- [MosekTools.jl](https://www.mosek.com/)
- [LinearAlgebra.jl](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- [Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/)
- [Random.jl](https://docs.julialang.org/en/v1/stdlib/Random/)
- [DataFrames.jl](https://dataframes.juliadata.org/stable/)
- [CSV.jl](https://csv.juliadata.org/stable/)
- [TypedPolynomials.jl](https://juliapackages.com/p/typedpolynomials)
- [MultivariateMoments.jl](https://juliapackages.com/p/multivariatemoments)
- [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/)
## Getting started
### completely positive matrix factorization bounds

``` julia
# Define the matrix of interest.
M =    [1.0       0.707107  0.0       0.0       0.447214
        0.707107  1.0       0.408248  0.0       0.0
        0.0       0.408248  1.0       0.288675  0.0
        0.0       0.0       0.288675  1.0       0.223607
        0.447214  0.0       0.0       0.223607  1.0]
# set the level of the hierarchy
lvl = 1 
# choose which constraitns to inclued (see the publication)
cons  = join(["G","dag","ddag","xx"][[1,2,4]],"")
# choose from the three main hierarchies (see the publication)
hier = ["id", "sp", "wsp"][2] 
# Define the model and run the optimization
ξₜᶜᵖⁱᵈ , _ = cp_model.get_ξₜᶜᵖ(M,lvl,hier*cons);
# Inspect the resutls in the terminal 
```
### attempt atom extraction
``` julia
mom_dir = assets_dir*"results\\lit\\t2\\moments\\"
mom_dir_list = [c for c in readdir(mom_dir) if contains(c,"id")]
mom_path = mom_dir*mom_dir_list[1]

  
df_id = matrix_IO.load_moments(mom_path)
include(src_dir*"extract_atoms.jl")

  
n_id, t_id, mom_vals_id = extract_atoms.proc_mom(df_id)
ext_atoms_id            = extract_atoms.get_atoms(n_id, t_id, mom_vals_id, 1e-4, true)
cents_id, weights_id    = extract_atoms.ext_centers_weights(ext_atoms_id)
A_ext_id                = extract_atoms.recon_mat(cents_id, weights_id)


```





