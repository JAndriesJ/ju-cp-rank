module batch_run
using LinearAlgebra
using Random
using CSV, DataFrames
Random.seed!(373)

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
           get_mom_mat_val_seq

    # Generates a batch of cp-matrices, stores them, their construciton, and rank 
    function gen_batch_sparse_cp_mats(n_range,num_ex,save_dir)
        for n in n_range, 
            p = rand()/2 +0.25
            for ex in num_ex
                ex = string(ex,pad=2)
                rp,R,M = cpm.get_random_sparse_cp_mats(n,p)
                r = LinearAlgebra.rank(M .+ 0.0)
                M_name = "M_n$(n)_ex$(ex)_r$(r)_rp$(rp).csv"
                cp_rank.save_mat(M, save_dir*M_name)
                R_name = "R_n$(n)_ex$(ex)_r$(r)_rp$(rp).csv"
                cp_rank.save_mat(R, save_dir*R_name)
            end
        end
    end
 
    # Loads all the matrices stored as .csv's from a directory as a dictionary
    function load_batch_sparse_cp_mats(load_dir)
        list_o_mats = readdir(load_dir)
        list_o_mat_names = [m[1:end-4] for m in list_o_mats]
        mats = [cp_rank.load_mat(load_dir*mat) for mat in list_o_mats]
        return Dict(zip(list_o_mat_names, mats))
    end

    function batch_comp_and_save(Mats_dict, t, save_dir)
        for M_name in sort([keys(Mats_dict)...])
            M = Mats_dict[M_name]
            n = size(M)[1]
            for cons in ["sGid","sG"]
                comp_save_path = save_dir*"ξₜᶜᵖ_$(M_name)_t$(t)_$cons.txt"
                ξₜᶜᵖ = run_and_save_get_ξₜᶜᵖ(M ,t ,cons, comp_save_path)
                mom_save_path = save_dir*"Moments_$(M_name)_t$(t)_$cons.csv"
                save_moments(ξₜᶜᵖ,n,t,mom_save_path)
            end
        end
    end

    function save_moments(ξₜᶜᵖ,n,t,save_path)
        df = DataFrame(moments_expo = ξₜᶜᵖ.obj_dict[:Lx].axes[1], 
                       moments_val = cp_model.rec_mom_mat(n,2*t,ξₜᶜᵖ))
                       CSV.write(save_path, df)  
    end

    load_moments(save_path) = make_mom_dic(CSV.read(save_path, DataFrame))
    make_mom_dic(df) = Dict(zip(clean_expo.(df.moments_expo), df.moments_val))
    clean_expo(expo) = parse.(Int,map(s -> strip(s,[',','[',']']),split(expo)))

    function make_summary(assets_dir)
        comp_dir = assets_dir*"comp\\"
        mat_dir = assets_dir*"mats\\"
        mom_dir = assets_dir*"moments\\"

        comp_names = readdir(assets_dir*"comp\\")
        file_loc = assets_dir*"summary.csv"
        create_summary_file(file_loc)

        Mats_dict = load_batch_sparse_cp_mats(mat_dir)

        for comp_name in comp_names
            mat_name = join(split(comp_name,'_')[2:end-2],'_')*".csv"
            mom_name = "Moments"*comp_name[12:end-4]*".csv"

            M = Mats_dict[mat_name[1:end-4]]
            mom_dict = load_moments(mom_dir*mom_name) 

            n, ex, r, rp ,t ,con = read_comp_file_name(comp_name) 
            sm, sms, smg, smgs = get_matrix_sizes(M,t)
            mom_rank = get_moment_ranks(n,t,mom_dict)
            get_matrix_sizes(M,t)
            try
                isP, isD, obj_v, sol_t = read_comp_file(comp_dir*comp_name)
                open(file_loc,"a") do io
                    write(io, "|$n|$ex|$r|$t|$con|$isP|$isD|$obj_v|$rp|$sol_t|$sm|$sms|$smg|$smgs|$mom_rank| \n")
                end
            catch
                open(file_loc,"a") do io
                    write(io, "|$n|$ex|$r|$t|$con|-|-|-|$rp|-|$sm|$sms|$smg|$smgs|$mom_rank|\n")
                end
            end
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

    function significant_sings(M) 
        u, s, v  = LinearAlgebra.svd(M)
        return  sum(cumsum(s) .<  0.99*cumsum(s)[end])
    end

    get_moment_ranks(n,t,mom_dict) = [significant_sings.(get_mom_mat_val_seq(n,t,mom_dict))]
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

    function create_summary_file(file_loc)
        touch(file_loc)
        open(file_loc,"a") do io
            write(io, "|n|ex|r|t|con|isP|isD|obj_v|rp|sol_t|sm|sms|smg|smgs|mom_ranks|\n")
        end
    end

    function read_comp_file_name(f_name)
        spl_name = split(f_name, '_')
        n =  parse(Int64,spl_name[3][2])
        ex = spl_name[4]
        r =  parse(Int64,spl_name[5][2:end])
        rp = parse(Int64,spl_name[6][3:end])
        t =  parse(Int64,spl_name[7][2])
        con = spl_name[8][1:end-4]
        return n, ex, r,rp ,t ,con
    end
    
    function read_comp_file(f_path)
        s = open(f_path) do file
            read(file, String)
        end
    
        nar = replace(s[end-95:end], '\n'=>' ')
        isPF = contains(nar,"Primal: FEASIBLE_POINT") 
        isDF = contains(nar,"Dual: FEASIBLE_POINT") 
    
        spl_kat = split(nar)
        sol_time = parse(Float64,spl_kat[3]) 
        obj_val = round(parse(Float64,spl_kat[9]), digits=5)  
    
        return isPF,isDF, obj_val,sol_time
    end
    
### -------------------------------------

    function clean_summary(assets_dir)
        f_name = assets_dir*"summary.csv"
        df = CSV.read(f_name,DataFrame,delim ="|")[:,2:end-1]
        df_obj_v = unstack(df, [:n,:ex,:r,:t,:rp,:sm,:sms,:smg,:smgs], :con, :obj_v, allowduplicates=true)
        rename!(df_obj_v , :sG => :sG_obj_v, :sGid => :sGid_obj_v)
        df_time  = unstack(df, [:n,:ex], :con, :sol_t, allowduplicates=true)
        rename!(df_time, :sG => :sG_time_s, :sGid => :sGid_time_s)
        df_mom_rank  = unstack(df, [:n,:ex], :con, :mom_ranks, allowduplicates=true)
        rename!(df_mom_rank, :sG => :sG_mom_ranks, :sGid => :sGid_mom_ranks)
        clean_df = innerjoin(df_obj_v, df_time, df_mom_rank, on = [:n,:ex])
        CSV.write(assets_dir*"clean_summary.csv", clean_df) 
    end



end