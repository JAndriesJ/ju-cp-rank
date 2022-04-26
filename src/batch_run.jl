module batch_run
using LinearAlgebra

using CSV, DataFrames
using FileIO


proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"cp_rank.jl")


using .cp_matrices ; const cpm = cp_matrices
using .moments ; const mom = moments
using .cp_model
using .cp_rank

    export gen_batch_sparse_cp_mats,
           load_batch_sparse_cp_mats,
           batch_comp_and_save,
           save_moments,
           load_moments,
           make_summary,
           make_summary_moments,
           clean_summary,
           get_mom_mat_val_seq,
           get_batch_mat_support

    function batch_comp_and_save(Mats_dict, t, G_act,save_dir)
        for M_name in sort([keys(Mats_dict)...])
            M = Mats_dict[M_name]
            for flavour ∈ ["ideal","gmpd","gmps"]
                try
                    comp_save_path = save_dir*"ξₜᶜᵖ_$(M_name)_$(flavour)_t$(t).txt"
                    cp_rank.get_ξₜᶜᵖ(M ,t, flavour, G_act, comp_save_path)
                catch
                end
                
            end

        end
    end
  ### -------------------------------------

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

    spl_kat = split(nar)
    sol_time = parse(Float64, spl_kat[findall(spl_kat .== "Time:")[1]+1])  
    obj_val  = round(parse(Float64,spl_kat[9]), digits=5)  

    return isPF, isDF, obj_val, sol_time
end

### -------------------------------------

function clean_summary(assets_dir)
    f_name = assets_dir*"summary.csv"
    df = CSV.read(f_name,DataFrame,delim ="|")[:,2:end-1]
    df_obj_v = unstack(df, [:name], :mod, :obj_v, allowduplicates=true)
    rename!(df_obj_v , :gmpd => :gmpd_obj_v, :gmps  => :gmps_obj_v, :ideal  => :ideal_obj_v)
    select!(df_obj_v,[:name, :ideal_obj_v, :gmpd_obj_v, :gmps_obj_v])

    df_time  = unstack(df, [:name], :mod, :sol_t, allowduplicates=true)
    rename!(df_time, :gmpd => :gmpd_time_s, :gmps => :gmps_time_s, :ideal => :ideal_time_s)
    select!(df_time,[:name, :ideal_time_s, :gmpd_time_s, :gmps_time_s])

    clean_df = innerjoin(df_obj_v, df_time, on = [:name])
    CSV.write(assets_dir*"clean_summary.csv", clean_df) 
end


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



    function get_batch_mat_support(load_dir,t=2)
        Mats_dict = batch_run.load_batch_sparse_cp_mats(load_dir)
        for mat_name in keys(Mats_dict)
            M = Mats_dict[mat_name]
            mat_supp = cp_rank.show_mat_support(M) ;
            mom_mat_supp = cp_rank.show_mat_support(moments.get_mom_mat_supp(t,M))
            ten_con_supp = cp_rank.show_mat_support(moments.get_ten_con_supp(t,M))   
    
            save(load_dir*"supp_"*mat_name*".png", mat_supp)
            save(load_dir*"mom_supp_t$t"*mat_name*".png", mom_mat_supp)
            save(load_dir*"ten_supp_t$t"*mat_name*".png", ten_con_supp)
        end
    end

    function get_matrix_sizes(M,t)
        n = size(M)[2]
        siz_Mₜ_spar =  length(mom.make_mon_expo(n,t,M))
        siz_Mₜ =length(mom.make_mon_expo(n,t))
        siz_MₜG_spar = length(mom.make_mon_expo(n,t-1,M))*n
        siz_MₜG = length(mom.make_mon_expo(n,t-1))*n
        return siz_Mₜ, siz_Mₜ_spar, siz_MₜG, siz_MₜG_spar
    end

### summary_moments-------------------------------------
    function make_summary_moments(save_dir)
        mom_names = readdir(save_dir)
        mom_ranks = []
        for mom in mom_names
            n,t = get_n_t(mom)
            mom_dict = load_moments(save_dir*mom) 
            mom_mat_rank_seq = get_moment_ranks(n,t,mom_dict)
            push!(mom_ranks,mom_mat_rank_seq)
        end
        return mom_ranks
    end

    get_moment_ranks(n,t,mom_dict) = [significant_sings.(get_mom_mat_val_seq(n,t,mom_dict))]
    function significant_sings(M) 
        u, s, v  = LinearAlgebra.svd(M)
        return  sum(cumsum(s) .<  0.99*cumsum(s)[end])
    end
    #get_moment_ranks(n,t,mom_dict) = [rank.(get_mom_mat_val_seq(n,t,mom_dict))]
    get_mom_mat_val_seq(n,t,mom_dict) = [α_to_Lxᵅ(mom_dict, mom_mat) for mom_mat ∈ get_mom_mat_seq(n,t)]
    get_mom_mat_seq(n,t) = [mom.make_mon_expo(n,(i,i)) for i ∈ 1:t]
    α_to_Lxᵅ(mom_dic::Dict{Vector{Int64}, Float64}, index_array) = map(α -> mom_dic[α], index_array)
    function get_n_t(mom_name) 
        n = parse(Int,split(mom_name,'_')[3][2])
        t = parse(Int,split(mom_name,'_')[end-1][2])
        return n,t
    end
    


### -------------------------------------   


    


end