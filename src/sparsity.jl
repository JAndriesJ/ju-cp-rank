module sparsity
using Test

export gen_sparsity,
       plot_matrix,
       plot_graph,
       gen_random_sparciy_mat,
       gen_random_sparciy_graph,
       clique_mem_num_vec 



using Random; Random.seed!(2017) # make sure this tutorial is reproducible
using LinearAlgebra

using Graphs # https://graphsjl-docs.readthedocs.io/en/latest/
# using GraphPlot, Plots, Fontconfig, Cairo, Colors


mutable struct sparsity_obj
    n::Int64
    V::Base.OneTo{Int64}
    E::Graphs.SimpleGraphs.SimpleEdgeIter{Graphs.SimpleGraphs.SimpleGraph{Int64}}
    adj_matrix::Matrix{Float64}
    supp_mat::Matrix{Int64}
    zero_ent::Vector{Tuple{Int64, Int64}}
    G::Graphs.SimpleGraphs.SimpleGraph{Int64}
    max_cliques::Vector{Vector{Int64}}
    max_clique_cover::Vector{Any} 
end


function gen_sparsity_obj(n,p=0.5)
    M = gen_random_sparciy_mat(n,p)
    G = Graphs.Graph(M)
    E = Graphs.edges(G)
    MC = Graphs.maximal_cliques(G)
    return sparsity_obj(n,
                        Graphs.vertices(G),
                        E,
                        Graphs.adjacency_matrix(G),
                        M,
                        get_zero_entries(M),
                        G,
                        MC,
                        get_max_clique_cover(E,MC))
end


### Plot functions


function plot_graph(G)
    E = Graphs.edges(G)
    mc = Graphs.maximal_cliques(G)
    mc = get_max_clique_cover(E,mc)
    cmnv = clique_mem_num_vec(E,mc)

    edgestrokec = distinguishable_colors(length(cmnv), colorant"blue")[cmnv] ;

    gplot(G, nodelabel=[Graphs.vertices(G)...],  
             edgestrokec=edgestrokec,
             layout=circular_layout)
end

### Matrix functions
bernoulli(p) =  rand() <= p ? 1 : 0
gen_random_sparciy_mat(n,p=0.5) = Int.(Symmetric(map(x -> bernoulli(p), ones(n,n))) + I(n) .> 0)
get_zero_entries(M) = [ (i,j) for i in 1:(size(M)[1]-1) for j in (i+1):size(M)[1] if M[i,j] == 0]


### Graph functions
gen_random_sparciy_graph(n,p=0.5) = Graphs.Graph(gen_random_sparciy_mat(n,p))

clique_mem_check(edge,clique) = edge.src ∈ clique && edge.dst ∈ clique
clique_mem(edge,clique_list) = map(c->clique_mem_check(edge,c),clique_list)
clique_mem_num(edge,clique_list) = findfirst(clique_mem(edge,clique_list))
clique_mem_num_vec(E,clique_list) = map(e->clique_mem_num(e,clique_list),E)

function get_max_clique_cover(E,MC)
    smc = MC[sortperm(length.(MC),rev=true)]
    clique_cover = []
    for e in E
        if any(clique_mem(e,clique_cover))
        else
            push!(clique_cover,smc[clique_mem_num(e,smc)])
        end
    end
    return clique_cover
end



function run_test()
    # Test 1
    g = graphfamous("karate")
    g = Graphs.smallgraph(:karate)
    gplot(g)


    # Test 2
    A = [0 1 1
         1 0 1
         1 1 0]
    G₂ = Graphs.Graph(A)
    gplot(G₂)


    # Test 4
    # E = Graphs.edges(G)
    # V = Graphs.vertices(G)
    # max_cli = Graphs.maximal_cliques(G)
    # Graphs.connected_components(G)
    # e = [E...][1]
    # c = max_cli[1]

    # clique_mem_check(e,c)
    # clique_mem(e,max_cli)
    # clique_mem_num(e,max_cli)
    # max_cli[12]
    # cmnv = clique_mem_num_vec(E,max_cli)


    M = [1 1 0 0 0 1
         1 1 1 1 0 1
         0 1 1 1 0 0
         0 1 1 1 1 1
         0 0 0 1 1 1
         1 1 0 1 1 1]
    
    G = Graphs.Graph(M)     
    Graphs.maximal_cliques(G)


    @testset "get_zero_entries" begin 
        up_tz = get_zero_entries(M)
        @test length(up_tz) == Int(sum(M .== 0)/2)
        @test [M[i,j] for (i,j) ∈ up_tz] == zeros(Int8,6)
    end
    
end


    
end