module extract_atoms

using DataFrames
using TypedPolynomials
using MultivariateMoments#, MultivariatePolynomials
using HomotopyContinuation

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")
include(proj_dir*"matrix_IO.jl")

export get_atoms

"""
Input
    mom_path: the file path of a .csv file containing the saved moments
    tol: the tolorance used in extracting atoms
    hom: Boolean, true for Homotopy continuation and false for Grobner basis. 
Output (either nothing or)
    cen : atom centers
    wei : atom weights 
    A   : The underlying matrix
"""
function get_atoms(mom_path::String, tol=1e-4, hom=true)
        # I am assuming a particular file structure here
        A = get_mat(mom_path) 

        df = matrix_IO.load_moments(mom_path)
        n, t, mom_vals = extract_atoms.proc_mom(df)
        ext_atoms = extract_atoms.get_atoms(n, t, mom_vals, tol, hom)
        try
                cen, wei = length(n) > 1 ? extract_atoms.ext_centers_weights(ext_atoms, A) : extract_atoms.ext_centers_weights(ext_atoms) 
                if length(n) > 1
                    A_rec = extract_atoms.recon_mat(cen,wei)
                else
                    A_rec = extract_atoms.recon_mat(cen)
                end
                
                return cen, wei, A, A_rec
        catch
                print("No atom could be extracted")
                return nothing, nothing, A, nothing
        end
end

function get_mat(mom_path)
    ex_type = split(mom_path,'\\')[end-3]
    mom_name = split(mom_path,'\\')[end]
    ex_name = join(split(mom_name,"_")[2:end-2],"_")
    ex_path = join(split(mom_path,'\\')[1:end-5],"\\")*"\\data\\"*ex_type*"\\"*ex_name*".csv"
    return  matrix_IO.load_mat(ex_path) 
end

"""(atoms,weights) -> recovered cp-matrix"""
recon_mat(c::Vector{Vector{Float64}})         = sum([c[i]*c[i]' for i in 1:length(c)])
recon_mat(c::Vector{Vector{Float64}}, w)      = sum([w[i]*c[i]*c[i]' for i in 1:length(c)])
recon_mat(c::Vector{Vector{Vector{Real}}}, w) = sum([ sum([w[j][i]*c[j][i]*c[j][i]' for i in 1:length(c[j])]) for j ∈ 1:length(c)])

"""Extracts the centers and weights of the atom objects"""
function ext_centers_weights(ext, A)
    n = size(A)[1] 
    mc = moments.get_maximal_cliques(A) 
    p = length(mc)
    cents   = [[[i ∈ mc[j] ? popfirst!(a.center) : 0  for i in 1:n] for a in ext[j].atoms] for j ∈ 1:p]
    weights = [[a.weight for a in ext[j].atoms] for j ∈ 1:p]
    return cents, weights
end
ext_centers_weights(extract) = [a.center for a in extract.atoms], [a.weight for a in extract.atoms]

"""takes the moment df and extracts: n, t, and moment values"""
proc_mom(df::DataFrame) = contains(string(df.mom_inds[1]),"]]") ? proc_sparse_mom(df) : proc_dense_mom(df)
function proc_dense_mom(df::DataFrame)
    n = length(df.mom_inds[1])
    t = div(maximum(maximum([m for m in df.mom_inds])),2) 
    mom_vals  = df.mom_vals 
    return  n, t, mom_vals
end
function proc_sparse_mom(df::DataFrame)
    mom_vals  = df.mom_vals 
    mom_inds = df.mom_inds
    nc = maximum([m[1] for m in mom_inds])
    n_s = [length([m[2] for m in mom_inds if m[1] == c][1]) for c in 1:nc]
    t = div(maximum([m[2][1] for m in mom_inds if m[1] == 1]),2) 
    val_s = [[ mom_vals[k] for k in 1:length(mom_vals) if mom_inds[k][1] == c] for c in 1:nc]
    return n_s, t, val_s
end

proc_mom(df::Matrix{Any}) = typeof(df[1]) == Vector{Any} ? proc_sparse_mom(df) : proc_dense_mom(df)
function proc_sparse_mom(df::Matrix{Any})
    m_i, m_v = get_inds_vals(df)
    p = length(m_v)
    nc  = maximum([m[1] for m in m_i])
    n_s = [length([m[2] for m in m_i if m[1] == c][1]) for c in 1:nc]
    t   = div(maximum([m[2][1] for m in m_i if m[1] == 1]),2) 
    val_s = [[ m_v[k] for k in 1:p if m_i[k][1] == c] for c in 1:nc]
    return n_s, t, val_s
end
function proc_dense_mom(df::Matrix{Any})
    m_i, m_v = get_inds_vals(df)
    n = length(mom_inds[1])
    t = div(maximum(maximum([m for m in m_i])),2) 
    return  n, t, m_v
end
get_inds_vals(df) = df[:,1],df[:,2]


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

