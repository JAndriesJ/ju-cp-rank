module cp_rank



proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_model.jl")

# using .cp_model  

export get_ξₜᶜᵖ   
## Run computations
function get_ξₜᶜᵖ(M,t,flavour,G_con=true) 
    if flavour=="id"
        mod = cp_model.modelξₜᶜᵖⁱᵈ(M,t,G_con=G_con)
    elseif flavour=="sp"
        mod = cp_model.modelξₜᶜᵖˢᵖ(M,t,G_con=G_con)
    elseif flavour=="wsp"
        mod = cp_model.modelξₜᶜᵖˢᵖ(M,t,G_con=G_con,isWeak=true)
    else
        error("Incorrect model specification")
    end
    return cp_model.computeξₜᶜᵖ(mod)
end  
function get_ξₜᶜᵖ(M ,t, flavour, save_path::String, G_con=true,)
    ξₜᶜᵖ, s = capture_solver_output(get_ξₜᶜᵖ, (M ,t, flavour, G_con))
    write_solver_output(s,save_path)
    ex_moments = cp_model.extract_moments(ξₜᶜᵖ)
    return ξₜᶜᵖ, ex_moments
end

function capture_solver_output(func,args)
    original_stdout = stdout;
    (rd, _) = redirect_stdout();
    #ξₜᶜᵖ = get_ξₜᶜᵖ(M ,t ,cons)
    ξₜᶜᵖ = func(args...)
    s = []
    for _ in 1:100
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

end


# using Plots, Fontconfig, Cairo
# using Graphs, GraphPlot

# include(proj_dir*"cp_matrices.jl")
# using .cp_matrices ; const cpm = cp_matrices

# get_spar_ξₜᶜᵖ,
# run_and_save_get_ξₜᶜᵖ,
# get_Bomze_cp_mats,
# get_random_cp_mat,
# get_random_sparse_cp_mat,
# gen_random_band_mat,
# show_mat_support,
# show_support_graph

# ### Visualizations
# function show_support_graph(M)
#     G = Graphs.Graph(M)
#     return gplot(G, nodelabel=[Graphs.vertices(G)...], layout=circular_layout)
# end
# show_mat_support(M,mag=1) = Gray.(repeat(M .> 0.0, inner=(mag,mag)))

# ## Generate matrices       

# get_random_cp_mat(n,r) = cpm.get_random_cp_mats(n,r)
# get_random_sparse_cp_mat(M::Matrix{Int64}) = cpm.get_random_sparse_cp_mats(M)
# get_random_sparse_cp_mat(n::Int,p::Float64=0.5) = cpm.get_random_sparse_cp_mats(cpm.gen_random_sparcity_mat(n,p))
# gen_random_band_mat(n,k) = cpm.gen_random_band_mat(n,k)