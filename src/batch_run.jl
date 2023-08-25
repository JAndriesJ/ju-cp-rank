module batch_run

using LinearAlgebra
using CSV, DataFrames
using FileIO
using Plots
using Plots.PlotMeasures

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"matrix_IO.jl")
include(proj_dir*"moments.jl")
include(proj_dir*"cp_model.jl")
include(proj_dir*"nn_model.jl")
include(proj_dir*"extract_atoms.jl")


using .moments ; const mom = moments
using .cp_model

export  batch_comp_and_save,
        make_summary,
        clean_summary,
        plot_bounds,
        plot_times,
        get_all_mom_mat_ranks


function get_ξₜ(args)
    if contains(args[3],"nn")
        ξₜ, ex_moments = nn_model.get_ξₜⁿⁿ(args...)
       
    elseif contains(args[3],"cp")
        ξₜ, ex_moments = cp_model.get_ξₜᶜᵖ(args...)
    else
        error("Must specify if the matrix is cp or just nn")
    end
    return ξₜ, ex_moments  
end

function get_ξₜ(args, save_path::String)
    ξₜ, mom, s = capture_solver_output(get_ξₜ, args) # args = (M ,t, flavour)
    write_solver_output(s, save_path)
    return ξₜ, mom
end
function capture_solver_output(func,args)
    original_stdout = stdout;
    (rd, _) = redirect_stdout();
    ξₜ, mom = func(args)
    s = []
    for rl in eachline(rd)
        push!(s,rl)
        if contains(rl,"Objective:")
            break
        end
    end
    redirect_stdout(original_stdout);
    return ξₜ, mom, s
end
function write_solver_output(s,save_path)
    touch(save_path)
    open(save_path, "w") do io
        for line in s
            write(io, line*"\n")
        end
    end;
end

function batch_comp_and_save(mats_dict, t, flavs, save_dir)
    mt = contains(flavs[1],"nn") ? "nn" : "cp"
    moments_dir = save_dir*"\\moments\\"
    !isdir(moments_dir) ? mkdir(moments_dir) : 0
    for n in sort([keys(mats_dict)...])
        M = mats_dict[n]
        for f ∈ flavs  
            comp_save_path = save_dir*"ξₜ$(mt)_$(n)_$(f)_t$(t).txt"
            mom_save_path  = moments_dir*"Mom_$(n)_$(f)_t$(t).csv"
            # try
                _ , ex_moments = get_ξₜ((M ,t, f), comp_save_path)
            matrix_IO.save_moments(ex_moments, mom_save_path)
            # catch
            # end
        end
    end
    return nothing
end

### Summarize data-------------------------------------

function make_summary(assets_dir)
    file_loc = assets_dir*"summary.csv"
    create_summary_file(file_loc)
    fs = readdir(assets_dir)
    tfns = fs[map(f -> contains(f,".txt"),fs)]
    for tf ∈ tfns
        n = join(split(tf,'_')[2:end-2],'_')
        m = split(tf,'_')[end-1]
        # try
            isP, isD, isUnkown, obj_v, sol_t = read_comp_file(assets_dir*tf)
            open(file_loc,"a") do io
                write(io, "$n,$m,$isP,$isD,$isUnkown,$obj_v,$sol_t \n")
            end
            println(tf*"-----------------success")
        # catch
        #     open(file_loc,"a") do io
        #         write(io, "$n,$m,,,, \n")
        #     end
        #     println(tf*"-----------------fail")
        # end
    end
    # ------------------Ranks of moment matrices------------------
    # mom_dir = assets_dir*"moments\\"
    # mom_ranks = [batch_run.get_mom_mat_ranks(mom_dir*m, 1e-4) for m in readdir(mom_dir)]
    # df_s = CSV.read(file_loc, DataFrame,delim =",")
    # df_s[!,"mom_ranks"] = mom_ranks
    # CSV.write(file_loc, df_s) 
    # ------------------Atoms (if they can be extracted) ------------------



end
function create_summary_file(file_loc)
    touch(file_loc)
    open(file_loc,"a") do io
        write(io, "name,mod,isP,isD,isUnkown,obj_v,sol_t\n")
    end
end
function read_comp_file(f_path)
    s = open(f_path) do file
        read(file, String)
    end
    ext_text  = replace(s[end-110:end], '\n'=>' ')
    isPF = contains(ext_text,"Primal: FEASIBLE_POINT") 
    isDF = contains(ext_text,"Dual: FEASIBLE_POINT") || contains(ext_text,"Dual: INFEASIBILITY_CERTIFICATE")
    isUnkown = false
    if contains(ext_text,"Primal: UNKNOWN_RESULT_STATUS") ||  contains(ext_text,"Dual: UNKNOWN_RESULT_STATUS")
        isPF = isDF = false
        isUnkown = true
    end

    spl  = split(ext_text)
    sol_time = round(parse(Float64, spl[findall(spl .== "Time:")[1]+1]), digits=2) 
    if (isPF && isDF && !isUnkown)
        obj_val  = round(parse(Float64,spl[end]), digits=2)
    elseif isUnkown
        obj_val  = "?" 
    else
        obj_val  = "*"
    end

    return isPF, isDF, isUnkown, obj_val, sol_time
end
#---------------------------------------------------------------------------------------
function clean_summary(results_dir)
    s_name = results_dir*"summary.csv"

    df_s = CSV.read(s_name, DataFrame,delim =",")
    
    df_obj_v = unstack(df_s, [:name], :mod, :obj_v, allowduplicates=true)
    df_time  = unstack(df_s, [:name], :mod, :sol_t, allowduplicates=true)

    df = innerjoin(df_obj_v, df_time, on = [:name],makeunique=true)
    
    CSV.write(results_dir*"clean_summary.csv", df) 
    return df
end

function clean_summary(results_dir, data_dir)
    d_name = data_dir*"mat_data.csv"
    s_name = results_dir*"summary.csv"

    df_s = CSV.read(s_name, DataFrame,delim =",")
    df_obj_v = unstack(df_s, [:name], :mod, :obj_v, allowduplicates=true)
    df_time  = unstack(df_s, [:name], :mod, :sol_t, allowduplicates=true)

    df_s = innerjoin(df_obj_v, df_time, on = [:name],makeunique=true)
    df_d = select!(CSV.read(d_name, DataFrame, delim =","), [:name, :ex ,:n, :nzd, :nc, :mc, :r, :ucpr])
    df = leftjoin(df_d, df_s , on = [:name])

    df.n = Int.(df.n)
    df.r = Int.(df.r)
    df.ucpr = Int.(df.ucpr)

    CSV.write(results_dir*"clean_summary.csv", df) 
    return df
end

### PLOTS -------------------------------------

function plot_bounds(ex_save_dir,t)
    df = CSV.read(ex_save_dir*"clean_clean_summary.csv", DataFrame,delim ="&")
    n = length(df.n)
    bounds =  names(df)[end-5:end-3]
    objs = scatter([1:n], Matrix(df[:,bounds]),  # , :ucpr
                                    xticks=(1:n, df.n),# df.ex .*"_n=".* string.(Int.(df.n))),
                                    yticks= Int(round(minimum(df[:,bounds[2]])-1)):Int(round(maximum(df[:,bounds[2]])+1)),
                                    yrange=[Int(round(minimum(df[:,bounds[3]])-1)),Int(round(maximum(df[:,bounds[2]])+1))],
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
    df = CSV.read(ex_save_dir*"clean_summary.csv", DataFrame,delim =",")
    n = length(df.n)
    times = names(df)[end-2:end]
    objs = scatter([1:n], Matrix(df[:,times]) .+ 1.0,  # , :ucpr
                                        xticks=(1:n, df.nzd),#df.ex .*"_nzd=".* string.(df.nzd)),# xlabel="Non-zero density for $(np)x$(np) matrices",
                                        xlabel="Size and nonzero density of the matrix",
                                        ylabel="time(s)",
                                        xrotation=90,
                                        left_margin = 20px,
                                        bottom_margin = 60px,
                                        label = ["id" "sp" "wsp"], # "e-rank_cp"
                                        legend=:topleft,
                                        fillalpha = 0.1,
                                        yscale=:log10,
                                        markershape = [:square :diamond :circle],
                                        markersizes = [4 4 3],
                                        markercolors = [:red :yellow :green],
                                        title = "Hierarchy times at t=$t",
                                        size = (1200, 600) )     #  :hline                                 
    savefig(objs, ex_save_dir*"times_t=$(t).png")    
                                   
end

function get_nn_ex_fig(results_dir)
        df_1 =  CSV.read(results_dir*"t1\\clean_summary.csv", DataFrame,delim =",")
        df_2 =  CSV.read(results_dir*"t2\\clean_summary.csv", DataFrame,delim =",")
        df_3 =  CSV.read(results_dir*"t3\\clean_summary.csv", DataFrame,delim =",")

        bounds_mat =  hcat(Matrix(df_1[:,["nn ddag id", "nn ddag sp"]]),
                        Matrix(df_2[:,["nn ddag id", "nn ddag sp"]]),
                        Matrix(df_3[:,["nn ddag id", "nn ddag sp"]]))
 
        disc_blocks =   2 .* Int.(bounds_mat' .>= 2) .+  Int.(bounds_mat' .> 3) .+  Int.(bounds_mat' .>= 4)            

        heatmap(disc_blocks,
                xrotation = 65,
                xticks=(0.5:5:100.5,0.0:0.05:1.0),
                yticks=(1.5:6.5, ["id t=1" "sp t=1" "id t=2" "sp t=2" "id t=3" "sp t=3"],),
                xaxis = "a",
                yaxis = "bound and level",
                # title = "ξₜⁿⁿˢᵖ and ξₜⁿⁿⁱᵈ bounds for t = 1,2,3 against a ∈ [0,1] ",
                c =cgrad([:yellow, :red, :purple], [2, 3, 4]),
                # c =cgrad([:white, :black], [0, 4]),
                size = (1000,200),
                bottom_margin = 10mm,
                left_margin = 10mm,
                legend=false,
                alpha = 0.7,
                gridalpha=1.0,
                grid=true,gridlinewidt=1.0)# ,gridstyle="solid"
  
end

#----------Moment rank checks-----------------   
function get_mom_mat_ranks(mom_path::String, ranktol =1e-4)
    df             = matrix_IO.load_moments(mom_path)
    n, t, mom_vals = extract_atoms.proc_mom(df)
    return get_mom_mat_ranks(n, t, mom_vals, ranktol)
end
get_mom_mat_ranks(n::Vector{Int64}, t, mom_vals,ranktol) = [get_mom_mat_ranks(n[j], t, mom_vals[j],ranktol) for j in 1:length(n)]
function get_mom_mat_ranks(n::Int64, t, mom_vals, ranktol =1e-4)
    r = []
    for j in 0:t
        mom_mat_thing  = extract_atoms.get_mom_mat(n, j, mom_vals)
        mom_mat        = [m for m in mom_mat_thing.Q]
        push!(r, lowrankchchek(mom_mat, ranktol))
    end
    return r
end
function lowrankchchek(M::AbstractMatrix, ranktol =1e-4)
    F = svd(M)
    nM = F.S[1] # norm of M
    tol = nM * ranktol
    r = something(findfirst(σ2 -> σ2 <= tol, F.S), 0)
    if iszero(r)
        cM = ranktol
        r = length(F.S)
    else
        cM = F.S[r] / nM
        r -= 1
    end
    r
end

end
