module cp_rank


using Plots, Fontconfig, Cairo
using Graphs, GraphPlot

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_matrices.jl")
include(proj_dir*"cp_model.jl")


using .cp_matrices ; const cpm = cp_matrices
using .cp_model  



export get_ξₜᶜᵖ,
       get_spar_ξₜᶜᵖ,
       run_and_save_get_ξₜᶜᵖ,
       get_Bomze_cp_mats,
       get_random_cp_mat,
       get_random_sparse_cp_mat,
       gen_random_band_mat,
       show_mat_support,
       show_support_graph
      
## Generate matrices       
get_Bomze_cp_mats() = cpm.get_Bomze_cp_mats()
get_random_cp_mat(n,r) = cpm.get_random_cp_mats(n,r)
get_random_sparse_cp_mat(M::Matrix{Int64}) = cpm.get_random_sparse_cp_mats(M)
get_random_sparse_cp_mat(n::Int,p::Float64=0.5) = cpm.get_random_sparse_cp_mats(cpm.gen_random_sparcity_mat(n,p))
gen_random_band_mat(n,k) = cpm.gen_random_band_mat(n,k)

## Run computations
function get_ξₜᶜᵖ(M,t,flavour,G_con) 
    if flavour=="ideal"
        mod = cp_model.modelξₜᶜᵖ(M,t,G_con=G_con)
    elseif flavour=="gmpd"
        mod = cp_model.GMP_dens_modelξₜᶜᵖ(M,t,G_con=G_con)
    elseif flavour=="gmps"
        mod = cp_model.GMP_spar_modelξₜᶜᵖ(M,t,G_con=G_con)
    else
        error("Incorrect model specification")
    end
    return cp_model.computeξₜᶜᵖ(mod)
end
    
function get_ξₜᶜᵖ(M ,t, flavour, G_con, save_path::String)
    ξₜᶜᵖ, s = capture_solver_output(get_ξₜᶜᵖ, (M ,t, flavour, G_con))
    write_solver_output(s,save_path)
    return ξₜᶜᵖ
end
function capture_solver_output(func,args)
    original_stdout = stdout;
    (rd, wr) = redirect_stdout();
    #ξₜᶜᵖ = get_ξₜᶜᵖ(M ,t ,cons)
    ξₜᶜᵖ = func(args...)
    s = []
    for i in 1:100
        rl = readline(rd)
        push!(s,rl)
        if contains(rl,"Objective:")
            break
        end
    end
    redirect_stdout(original_stdout);
    return ξₜᶜᵖ, s
end
function write_solver_output(s,save_path)
    touch(save_path)
    open(save_path, "w") do io
        for line in s
            write(io, line*"\n")
        end
    end;
end




### Visualizations
function show_support_graph(M)
    G = Graphs.Graph(M)
    return gplot(G, nodelabel=[Graphs.vertices(G)...], layout=circular_layout)
end
show_mat_support(M,mag=1) = Gray.(repeat(M .> 0.0, inner=(mag,mag)))


function run_tests()
    cpm.run_tests()
    mom.run_tests()
    con.run_tests()
    cp_model.run_tests()  
    compute.run_tests() 
end

end
