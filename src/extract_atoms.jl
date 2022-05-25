module extract_atoms

using TypedPolynomials
using MultivariateMoments#, MultivariatePolynomials
using HomotopyContinuation

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"batch_run.jl")
using .moments 
using .batch_run 

export proc_mom,
       get_mom_mat,
       get_atoms,
       ext_centers_weights,
       recon_mat

"""(atoms,weights) -> recovered cp-matrix"""
recon_mat(c::Vector{Vector{Float64}}) = sum([c[i]*c[i]' for i in 1:length(c)])
recon_mat(c::Vector{Vector{Float64}}, w) = sum([w[i]*c[i]*c[i]' for i in 1:length(c)])
recon_mat(c::Vector{Vector{Vector{Real}}}, w) = sum([ sum([w[j][i]*c[j][i]*c[j][i]' for i in 1:length(c[j])]) for j ∈ 1:length(c)])

"""Extracts the centers and weights of the atom objects"""
function ext_centers_weights(extract, A)
    n = size(A)[1] ; mc = moments.get_maximal_cliques(A) 
    cents_sp = [[[i ∈ mc[j] ? popfirst!(a.center) : 0  for i in 1:n] for a in extract[j].atoms] for j ∈ 1:length(mc)]
    weights_sp = [[a.weight for a in extract[j].atoms] for j ∈ 1:length(mc)]
    return cents_sp, weights_sp
end
function ext_centers_weights(extract)
    cents = [a.center for a in extract.atoms]
    weights = [a.weight for a in extract.atoms]
    return cents, weights
end

"""takes the moment df and extracts: n, t, and moment values"""
function proc_mom(df)
    if contains(df.mom_inds[1],"]]")
        return proc_sparse_mom(df)
    else
        return proc_dense_mom(df)
    end
end
function proc_dense_mom(df)
    mom_vals  = df.mom_vals 
    mom_inds = [eval(Meta.parse(ind)) for ind in df.mom_inds]
    n = length(mom_inds[1])
    t = div(maximum(maximum([m for m in mom_inds])),2) 
    return  n, t, mom_vals
end
function proc_sparse_mom(df)
    mom_vals  = df.mom_vals 
    mom_inds = [eval(Meta.parse(ind)) for ind in df.mom_inds]
    nc = maximum([m[1] for m in mom_inds])
    n_s = [length([m[2] for m in mom_inds if m[1] == c][1]) for c in 1:nc]
    t = div(maximum([m[2][1] for m in mom_inds if m[1] == 1]),2) 
    val_s = [[ mom_vals[k] for k in 1:length(mom_vals) if mom_inds[k][1] == c] for c in 1:nc]
    return n_s, t, val_s
end

###
"""Builds moment matrix and extracts """
get_atoms(n::Vector{Int64}, t, mom, tol=1e-4, hom=true) = [get_atoms(n[j],t,mom[j], tol,hom) for j in 1:length(n) ]
function get_atoms(n::Int64, t, mom, tol=1e-4, hom=true)
    M = get_mom_mat(n, t, mom)
    if hom 
        solver = HomotopyContinuation.SemialgebraicSetsHCSolver(; compile = true)
        return extractatoms(M, tol, solver)
    else
        return extractatoms(M, tol)
    end
end

get_mom_mat(n::Vector{Int64}, t::Int64, mom) = [get_mom_mat(n[j], t, mom[j]) for j in 1:length(n)]
function get_mom_mat(n::Int64, t::Int64, mom)
    x, μ = get_psudomeasure(n, t, mom)
    B = map(a -> prod(x.^a), moments.make_mon_expo(n, t))
    return MultivariateMoments.moment_matrix(μ, B)
end

function get_psudomeasure(n,t,mom_vals)
    x, monomials = get_psudome_set_up(n,t)
    return x, MultivariateMoments.measure(mom_vals[1:length(monomials)] , monomials)
end
function get_psudome_set_up(n,t)
    @HomotopyContinuation.polyvar x[1:n]       
    monomials_expos = moments.make_mon_expo(n, 2t)
    monomials = map(a -> prod(x.^a), monomials_expos)
    return x, monomials
end



end

# function get_psudomeasure(n,t,load_path::String)
#     mom_vals = batch_run.load_moments(load_path)
#     x, monomials = get_psudome_set_up(n,t)
#     return x, MultivariateMoments.measure(mom_vals , monomials)
# end