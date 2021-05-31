#"""save computations to a .txt file """
module Batch_proc
using JuMP


proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cpMatrices.jl")
include(proj_dir*"cpModel.jl")
include(proj_dir*"Compute.jl")

include("C:\\Users\\andries\\all-my-codes\\ju-text-manipulation-tools\\read_n_proc_tools.jl")
using .read_n_proc_tools

using .cpMatrices
using .cpModel
using .Compute

export batchModelξ₂ᶜᵖ,
       batchCompξᶜᵖ

function gen_con_names()
    DagXX = vec(string.(["","Dag"],["" "XX"]))
    con_list = sort(vec(string.(DagXX ,["" "wG" "sG"])))
    return con_list
end

function batchModelξ₂ᶜᵖ(t::Int,save_dir::String)
    counter = 1
    con_list = gen_con_names()
    for MatName in ["M6", "M7", "M7t", "M8t", "M9t","M11t"]
        A = gen_cp_mats(MatName)
        for con in []
            println("############################ Matrix: $MatName Constraints: $con #############################  ")

            Mat_con_Name = MatName*"_"*con
            model = Modelξₜᶜᵖ(A,t,con)
            write_to_file(model, save_dir*"$Mat_con_Name.dat-s", format=MOI.FileFormats.FORMAT_SDPA)
        end
    end
end

function batchCompξᶜᵖ(load_dir)
    dats_files = collect_files_in_dir(load_dir, ".dat-s")
    cdir = split(load_dir,"\\")[end-1]

    filename = load_dir*cdir*".md"
    touch(filename)
    for dats_file in dats_files
        dats_file_name = split(dats_file,".")[end-1]
        con = split(dats_file_name,cdir)[end]

        if con !=  "DagXXsG" # "-xx-dag-G"
            continue
        end

        if con == ""
            con = "-"
        end

        cp_model = read_dat_s_model(dats_file)
        cp_model_opt = Computeξₜᶜᵖ(cp_model)
        pdo = extract_model_stats(cp_model_opt)

        open(filename, "a") do f
            write(f, "$con, $pdo  \n")
        end
    end
end

end
