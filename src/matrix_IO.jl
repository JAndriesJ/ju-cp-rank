module matrix_IO
using CSV, DataFrames

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"moments.jl")

### Matrices
save_mat(M, save_path::String) = CSV.write(save_path, DataFrame(M, :auto))
load_mat(load_path::String) = Matrix(CSV.read(load_path, DataFrame))
"""Loads all the matrices stored as .csv's from a directory as a dictionary"""
function load_mats(load_dir::String)
    if load_dir[end-3:end] == ".csv"
      return load_mat(load_dir) 
    else
      list_o_mats = [s for s in readdir(load_dir) if (contains(s,".csv") && contains(s,"ex"))]
      list_o_mat_names = [m[1:end-4] for m in list_o_mats]
      mats = [load_mat(load_dir*mat) for mat in list_o_mats]
      return sort(Dict(zip(list_o_mat_names, mats)))
    end
  end

### Moments
function save_moments(ext_mom, save_path::String)  
    df = DataFrame( mom_inds = ext_mom[:,1],
                    mom_vals = ext_mom[:,2]) 
    CSV.write(save_path, df)
end

function load_moments(load_path::String)
    df = CSV.read(load_path, DataFrame)
    df.mom_inds = eval.(Meta.parse.(df.mom_inds))
    return df
end
### String manipulation
function split_name(name::String) 
  sname= split(name,'_')
  n    = parse(Int64, sname[2][2:end])
  zd   = parse(Float64,sname[3][4:end])
  r    = parse(Int64, sname[4][2:end])
  ucpr = parse(Int64, sname[5][5:end])
  return [n, zd, r, ucpr]
end

function get_mat_data(data_dir::String)
  mats = matrix_IO.load_mats(data_dir)
  K = [keys(mats)...]
  mcs = map(k -> moments.get_maximal_cliques(mats[k]), K) 
  max_cliq = [maximum(length.(mc)) for mc in mcs ]
  sn = cat(map(name->matrix_IO.split_name(name), K)..., dims=2)'
  df = DataFrame(     name = K,
                      ex   = map(name -> name[1:4],K), 
                      n    = sn[:,1],
                      nzd  = sn[:,2],
                      r    = sn[:,3],
                      ucpr = sn[:,4],
                      nc   = length.(mcs),
                      mc   = max_cliq)
 
  CSV.write(data_dir*"mat_data.csv", df)
end

end

