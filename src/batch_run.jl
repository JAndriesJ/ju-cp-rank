module batch_run

using LinearAlgebra
using CSV, DataFrames
using FileIO
using Plots
using Plots.PlotMeasures

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"cp_rank.jl")

using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .cp_model
using .cp_rank

export  get_mat_data,
        batch_comp_and_save,
        make_summary,
        clean_summary,
        clean_clean_summary,
        plot_bounds,
        plot_times,
        save_moments,
        load_moments

function batch_comp_and_save(Mats_dict, t, save_dir, flavs,G_act=true)
    for M_name in sort([keys(Mats_dict)...])
        M = Mats_dict[M_name]
        for fla ∈ flavs  # 
            comp_save_path = save_dir*"ξₜᶜᵖ_$(M_name)_$(fla)_t$(t).txt"
            _ , ex_moments = cp_rank.get_ξₜᶜᵖ(M ,t, fla, comp_save_path,G_act)
            mom_save_path = save_dir*"Mom_$(M_name)_$(fla)_t$(t).csv"
            save_moments(ex_moments, mom_save_path)
        end
    end
end

### Summarize data-------------------------------------

function make_summary(assets_dir)
    file_loc = assets_dir*"summary.csv"
    create_summary_file(file_loc)

    fs = readdir(assets_dir)
    text_file_names = fs[map(f -> contains(f,".txt"),fs)]
    for comp_name in text_file_names
        mat_name = join(split(comp_name,'_')[2:end-2],'_')
        println(comp_name)
        mod = split(comp_name,'_')[end-1]
        try
            isP, isD, obj_v, sol_t = read_comp_file(assets_dir*comp_name)
            open(file_loc,"a") do io
                write(io, "|$mat_name|$mod|$isP|$isD|$obj_v|$sol_t| \n")
            end
        catch
            open(file_loc,"a") do io
                write(io, "|$mat_name|$mod|0|0|0|0| \n")
            end
        end
    end
end
function create_summary_file(file_loc)
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "|name|mod|isP|isD|obj_v|sol_t|\n")
    end
end
function read_comp_file(f_path)
    s = open(f_path) do file
        read(file, String)
    end

    nar = replace(s[end-95:end], '\n'=>' ')
    isPF = contains(nar,"Primal: FEASIBLE_POINT") 
    isDF = contains(nar,"Dual: FEASIBLE_POINT") 

    spl_kat  = split(nar)
    sol_time = parse(Float64, spl_kat[findall(spl_kat .== "Time:")[1]+1])  
    obj_val  = round(parse(Float64,spl_kat[9]), digits=5)  

    return isPF, isDF, obj_val, sol_time
end

function clean_summary(assets_dir)
    f_name = assets_dir*"summary.csv"
    df = CSV.read(f_name,DataFrame,delim ="|")[:,2:end-1]
    df_obj_v = unstack(df, [:name], :mod, :obj_v, allowduplicates=true)
    rename!(df_obj_v , :id => :id_obj_v, :sp => :sp_obj_v, :wsp  => :wsp_obj_v )
    select!(df_obj_v,[:name, :id_obj_v, :sp_obj_v, :wsp_obj_v]) 

    df_time  = unstack(df, [:name], :mod, :sol_t, allowduplicates=true)
    rename!(df_time, :id => :id_time_s, :sp => :sp_time_s, :wsp => :wsp_time_s)
    select!(df_time,[:name, :id_time_s, :sp_time_s, :wsp_time_s])

    df = innerjoin(df_obj_v, df_time, on = [:name])
    df.id_obj_v  = round.(df.id_obj_v  ,digits=2)
    df.sp_obj_v  = round.(df.sp_obj_v  ,digits=2) 
    df.wsp_obj_v = round.(df.wsp_obj_v ,digits=2)

    CSV.write(assets_dir*"clean_summary.csv", df) 
end

function clean_clean_summary(assets_dir, data_file)
    f_name = assets_dir*"clean_summary.csv"
    df_1 = CSV.read(data_file, DataFrame, delim ="&")
    df_2 = CSV.read(f_name, DataFrame, delim =",")

    df = innerjoin(df_1, df_2 , on = [:name])
    select!(df, [:ex ,:n, :nzd, :nc, :mc, :r, :id_obj_v, :sp_obj_v, :wsp_obj_v, :ucpr , :id_time_s, :sp_time_s, :wsp_time_s])

    df.n = Int.(df.n)
    df.r = Int.(df.r)
    df.ucpr = Int.(df.ucpr)

    df.id_obj_v = round.(df.id_obj_v,digits=2) 
    df.sp_obj_v = round.(df.sp_obj_v,digits=2) 
    df.wsp_obj_v = round.(df.wsp_obj_v,digits=2) 

    df.id_time_s = round.(df.id_time_s,digits=2) 
    df.sp_time_s = round.(df.sp_time_s,digits=2) 
    df.wsp_time_s = round.(df.wsp_time_s,digits=2)  

    CSV.write(assets_dir*"clean_clean_summary.csv", df,delim="&") 
end
function split_name(name) 
    dar = split(name,'_')
    ex   = dar[1]
    n    = parse(Int64,dar[2][2:end])
    zd   = parse(Float64,dar[3][4:end])
    r    = parse(Int64,dar[4][2])
    ucpr = parse(Int64,dar[5][5:end])
    return [n, zd, r, ucpr]
end

function get_mat_data(data_dir)
    mats = cp_matrices.load_mats(data_dir)
    K = [keys(mats)...]
    mcs = map(k -> moments.get_maximal_cliques(mats[k]), K)
    num_cliq = length.(mcs) 
    max_cliq = [maximum(length.(mc)) for mc in mcs ]
    dara = cat(map(name->split_name(name), [keys(mats)...])..., dims=2)'
    df = DataFrame(     name = K,
                        ex   = map(name -> name[1:4],K), 
                        n    = dara[:,1],
                        nzd  = dara[:,2],
                        r    = dara[:,3],
                        ucpr = dara[:,4],
                        nc   = num_cliq,
                        mc   = max_cliq)
   
    CSV.write(data_dir*"mat_data.csv", df, delim="&")
end
get_z_dense(M) = sum(M .== 0) / ((size(M)[1])*(size(M)[1]-1))

### PLOTS -------------------------------------

function plot_bounds(ex_save_dir,t)
    df = CSV.read(ex_save_dir*"clean_clean_summary.csv", DataFrame,delim ="&")
    n = length(df.n)
    np = df.n[1]
    objs = scatter([1:n], Matrix(df[:,[:id_obj_v, :sp_obj_v, :wsp_obj_v]]),  # , :ucpr
                                    xticks=(1:n, df.n),# df.ex .*"_n=".* string.(Int.(df.n))),
                                    yticks= Int(round(minimum(df.sp_obj_v)-1)):Int(round(maximum(df.sp_obj_v)+1)),
                                    yrange=[Int(round(minimum(df.wsp_obj_v)-1)),Int(round(maximum(df.sp_obj_v)+1))],
                                    # xlabel="non-zero density for $(np)x$(np) matrices",
                                    xlabel="Size of the matrix",
                                    ylabel="bound",
                                    xrotation=90,
                                    left_margin = 20px,
                                    bottom_margin = 60px,
                                    label = ["ξₜᶜᵖ⁻ⁱᵈ" "ξₜᶜᵖ⁻ˢᵖ" "ξₜᶜᵖ⁻ʷˢᵖ"], # "e-rank_cp"
                                    legend=:topleft,
                                    fillalpha = 0.1,
                                    markershape = [:square :diamond :circle],
                                    markersizes = [6 6 4],
                                    markercolors = [:red :yellow :green],
                                    title = "Hierarchy bounds at t=$t",
                                    size = (900, 600) )     #  :hline  
    savefig(objs, ex_save_dir*"bound_t=$(t).png")
end

function plot_times(ex_save_dir,t)
    df = CSV.read(ex_save_dir*"clean_clean_summary.csv", DataFrame,delim ="&")
    n = length(df.n)
    np = df.n[1]
    objs = scatter([1:n], Matrix(df[:,[:id_time_s, :sp_time_s, :wsp_time_s]]) .+ 1.0,  # , :ucpr
                                        xticks=(1:n, df.n),#df.ex .*"_nzd=".* string.(df.nzd)),# xlabel="Non-zero density for $(np)x$(np) matrices",
                                        xlabel="Size of the matrix",
                                        ylabel="time(s)",
                                        xrotation=90,
                                        left_margin = 20px,
                                        bottom_margin = 60px,
                                        label = ["ξₜᶜᵖ⁻ⁱᵈ" "ξₜᶜᵖ⁻ˢᵖ" "ξₜᶜᵖ⁻ʷˢᵖ"], # "e-rank_cp"
                                        legend=:topleft,
                                        fillalpha = 0.1,
                                        yscale=:log10,
                                        markershape = [:square :diamond :circle],
                                        markersizes = [6 6 4],
                                        markercolors = [:red :yellow :green],
                                        title = "Hierarchy times at t=$t",
                                        size = (900, 600) )     #  :hline 
    savefig(objs, ex_save_dir*"times_t=$(t).png")                                    
end


### Moments
function save_moments(ext_mom, save_path)  
    df = DataFrame( mom_inds = ext_mom[:,1],
                    mom_vals = ext_mom[:,2]) 
    CSV.write(save_path, df)
end

"""moment matrix M₂ₜ(y)"""
function load_moments(n::Int, t::Tuple{Int,Int}, load_path::String) 
    mom_vals = load_moments(load_path).mom_vals
    mom_dict = make_mom_dict(n::Int,t[1], mom_vals)
    return map(a -> mom_dict[a], make_mon_expo(n, (t[1],t[2])))
end
"""Dictionary of monomial exponent β ∈ ℕⁿ₂ₜ keys and values yᵦ"""
function make_mom_dict(n::Int,t::Int, mon_vals)
    # @assert length(mon_vals) ....
    mon_expo = make_mon_expo(n,2t)
    Dict(zip(mon_expo, mon_vals))
end
load_moments(load_path) = CSV.read(load_path, DataFrame)


#---------------------------   
       


end

### Exploratory -------------------------------------

# function get_batch_mat_support(load_dir,t=2)
#     Mats_dict = batch_run.load_batch_sparse_cp_mats(load_dir)
#     for mat_name in keys(Mats_dict)
#         M = Mats_dict[mat_name]
#         mat_supp = cp_rank.show_mat_support(M) ;
#         mom_mat_supp = cp_rank.show_mat_support(moments.get_mom_mat_supp(t,M))
#         ten_con_supp = cp_rank.show_mat_support(moments.get_ten_con_supp(t,M))   

#         save(load_dir*"supp_"*mat_name*".png", mat_supp)
#         save(load_dir*"mom_supp_t$t"*mat_name*".png", mom_mat_supp)
#         save(load_dir*"ten_supp_t$t"*mat_name*".png", ten_con_supp)
#     end
# end

# function get_matrix_sizes(M,t)
#     n = size(M)[2]
#     siz_Mₜ_spar =  length(mom.make_mon_expo(n,t,M))
#     siz_Mₜ =length(mom.make_mon_expo(n,t))
#     siz_MₜG_spar = length(mom.make_mon_expo(n,t-1,M))*n
#     siz_MₜG = length(mom.make_mon_expo(n,t-1))*n
#     return siz_Mₜ, siz_Mₜ_spar, siz_MₜG, siz_MₜG_spar
# end

# function plot_summary(save_dir)
#     df = CSV.read(save_dir*"clean_clean_summary.csv",DataFrame,delim =",")
#     dara = cat(map(name->split_name(name), df.name )..., dims=2)'
#     df_n = DataFrame(name = df.name,
#                         n    = dara[:,1],
#                         p    = dara[:,2],
#                         r    = dara[:,3],
#                         ucpr = dara[:,4])
#     plot_df = select!(innerjoin(df_n, df, on = [:name]),[:n, :r, :id_obj_v,  :sp_obj_v,  :wsp_obj_v, :ucpr])
#     n = size(plot_df)[1]
#     objs = scatter([1:n], Matrix(plot_df[:,[:r, :id_obj_v, :sp_obj_v, :wsp_obj_v]]),  # , :ucpr
#                                 xticks=(1:n,plot_df.n),
#                                 yticks=(5:15,5:15),
#                                 xlabel="size",
#                                 ylabel="bound",
#                                 xrotation=90,
#                                 # bottom_margin = 130px,
#                                 label = ["rank" "ξₜᶜᵖ⁻ⁱᵈ" "ξₜᶜᵖ⁻ˢᵖ" "ξₜᶜᵖ⁻ʷˢᵖ"], # "e-rank_cp"
#                                 legend=:topleft,
#                                 fillalpha = 0.1,
#                                 markershape = [:hline :xcross :cross :rtriangle],
#                                 size = (900, 600) )     #  :hline  
#     return objs
# end   
# function split_name(name)
#     spname = split(name,'_')
#     n = parse(Int, spname[2][2:end])
#     zd = parse(Float64, spname[3][3:end])
#     r = parse(Int, spname[4][2:end])
#     ucpr = parse(Int, spname[5][5:end])
#     return [n,zd,r,ucpr]
# end


### -------------------------------------
    # mom_save_path = save_dir*"Moments_$(M_name)_t$(t)_$cons.csv"
    # save_moments(ξₜᶜᵖ,n,t,mom_save_path)

# function read_comp_file_name(f_name)
#     spl_name = split(f_name, '_')
#     n   =  parse(Int64, spl_name[3][2])
#     ex  =  spl_name[4]
#     r   =  parse(Int64, spl_name[5][2:end])
#     rp  =  parse(Int64, spl_name[6][3:end])
#     t   =  parse(Int64, spl_name[7][2])
#     con =  spl_name[8][1:end-4]
#     return n, ex, r, rp ,t ,con
# end

# n, ex, r, rp ,t ,con = read_comp_file_name(comp_name) 
# mat_dir  = assets_dir*"mats\\"
    # mom_dir  = assets_dir*"moments\\"
    # Mats_dict = load_batch_sparse_cp_mats(mat_dir)
# sm, sms, smg, smgs = get_matrix_sizes(Mats_dict[mat_name[1:end-4]], t)
#mom_name = "Moments"*comp_name[12:end-4]*".csv"
#mom_rank = get_moment_ranks(n,t,load_moments(mom_dir*mom_name))    
### -------------------------------------

### summary_moments-------------------------------------
# function make_summary_moments(save_dir)
#     mom_names = readdir(save_dir)
#     mom_ranks = []
#     for mom in mom_names
#         n,t = get_n_t(mom)
#         mom_dict = load_moments(save_dir*mom) 
#         mom_mat_rank_seq = get_moment_ranks(n,t,mom_dict)
#         push!(mom_ranks,mom_mat_rank_seq)
#     end
#     return mom_ranks
# end

# get_moment_ranks(n,t,mom_dict) = [significant_sings.(get_mom_mat_val_seq(n,t,mom_dict))]
# function significant_sings(M) 
#     u, s, v  = LinearAlgebra.svd(M)
#     return  sum(cumsum(s) .<  0.99*cumsum(s)[end])
# end
# #get_moment_ranks(n,t,mom_dict) = [rank.(get_mom_mat_val_seq(n,t,mom_dict))]
# get_mom_mat_val_seq(n,t,mom_dict) = [α_to_Lxᵅ(mom_dict, mom_mat) for mom_mat ∈ get_mom_mat_seq(n,t)]
# get_mom_mat_seq(n,t) = [mom.make_mon_expo(n,(i,i)) for i ∈ 1:t]
# α_to_Lxᵅ(mom_dic::Dict{Vector{Int64}, Float64}, index_array) = map(α -> mom_dic[α], index_array)
# function get_n_t(mom_name) 
#     n = parse(Int,split(mom_name,'_')[3][2])
#     t = parse(Int,split(mom_name,'_')[end-1][2])
#     return n,t
# end

# ### ----------