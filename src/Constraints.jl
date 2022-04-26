module constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  α_to_Lxᵅ,
        make_dag_con,
        make_loc_con,
        make_xx_con,
        make_G_con,
        make_ideal_con
       
"""
input: A(Data matrix),t(integer),Lx(JuMP variable)
output: dictionary: keys:
L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)
"""
function make_dag_con(A,t,Lx)
    n = size(A)[1]
    dag_con = Dict()

    mom₂ₜ     = mom.make_mon_expo(n,2*t)
    dag_con[(0,0)] = [ 1*Lx[α ] for α in mom₂ₜ] # This is the non negativity of the moments i.e. products of x's

    deg_g = 2
    mom₂ₜ₋₂    = mom.make_mon_expo(n,2*t-deg_g)
    @assert lastindex(mom₂ₜ₋₂) > 1 "Contraints do not exist for t =1!"
    for k in 1:n
        eₖ = mom.eᵢ(n,k)
        sqrMₖₖ = sqrt(A[k,k])
        dag_con[(k,k)] = [ sqrMₖₖ*Lx[eₖ + α] - Lx[2*eₖ + α] for α in mom₂ₜ₋₂] # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        for h in (k+1):n
            eₕ = mom.eᵢ(n,h)
            dag_con[(k,h)] = [ A[k,h]*Lx[α] - Lx[eₖ + eₕ + α]  for α in mom₂ₜ₋₂] # Dagger constraints: L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        end
    end
    return dag_con
end

"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (i,j) ∈ [n]²
                    vals:
comment: L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ)
(M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function make_loc_con(A,t,Lx)
    n = size(A)[1]
    nze = mom.get_nonzero_entries(A)
    momₜ₋₁  = mom.make_mon_expo(t-1, A)
    # if_diag(k) = [sqrt(A[k,k])*Lx[α+β+mom.eᵢ(n,k)] - Lx[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    # if_off_diag(k,h) = [A[k,h]*Lx[α+β] - Lx[α+β+mom.eᵢⱼ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    momₜ₋₁g(j)  = mom.make_mon_expo(t-1,A,j)
    if_diag(k) = [sqrt(A[k,k])*Lx[α+β+mom.eᵢ(n,k)] - Lx[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁g(k), β ∈ momₜ₋₁g(k)]
    if_off_diag(k,h) = [A[k,h]*Lx[α+β] - Lx[α+β+mom.eᵢ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    return [e[1]==e[2] ? if_diag(e[1]) : if_off_diag(e[1],e[2]) for e in nze]
end

"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (h,k) ∈ [n]², h ≠ k
                    values:
comment: L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(xₖxₕL) )_αβ =   x^(α + β + eₖ + eₕ) 
"""
function make_xx_con(A,t,Lx)
    n       = size(A)[1]
    momₜ₋₁   = mom.make_mon_expo(n,t-1)
    eᵢⱼs    = map( e -> mom.eᵢ(n,e[1],e[2]), mom.get_nonzero_entries(A)) 
    return map(eᵢⱼ -> [Lx[α+β+eᵢⱼ] for α ∈ momₜ₋₁, β ∈ momₜ₋₁ ], eᵢⱼs)
end

"""M(G ⊗ L) ⪰ 0 constraints
input: A(data matrix),t(Integer),Lx(JuMP variable)
Assumption: G = A-[x]₌₁[x]₌₁ᵀ
output: A⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
"""
function make_G_con(A,t,Lx)
    n = size(A)[1]
    LMBexp₁    = mom.make_mon_expo(n,(1,1),isle=false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexpₜ₋₁   = mom.make_mon_expo((t-1,t-1),A)     #exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁      = α_to_Lxᵅ(Lx,LMBexpₜ₋₁)    #L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp₁ₜ₋₁  = mom.expo_kron(LMBexp₁,LMBexpₜ₋₁)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB₁ₜ₋₁     = α_to_Lxᵅ(Lx,LMBexp₁ₜ₋₁)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB₁ₜ₋₁             # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    return G_con
end

"""M(xᵢxⱼL) = 0  for all {i,j} s.t. Mᵢⱼ = 0"""
function make_ideal_con(A,t,Lx) 
    n = size(A)[1] 
    momₜ₋₁ = mom.make_mon_expo(t-1,A)
    eᵢⱼs = map( e -> mom.eᵢ(n,e[1],e[2]), mom.get_zero_entries(A))
    return map(eᵢⱼ -> [Lx[α+β+eᵢⱼ] for α ∈ momₜ₋₁, β ∈ momₜ₋₁ ], eᵢⱼs)
end

############################################### Utility #########################################
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Lx, index_array) = map(α -> Lx[α], index_array)


function run_tests()
    @testset "sparse" begin
        # @test var_inds = get_var_inds(mc,t) 
        # @test sort(unique([ b  for (a,b) ∈ var_inds])) ==  sort(mom.make_mon_expo(length(mc[k]),2*t))
        # @test length(spar_pos_cons) == length(cat(get_maximal_cliques(M)...,dims=1))
        # @test length(spar_loc_cons) == length(cat([nz(m) for m in get_maximal_cliques(M)]...,dims=1))
    end
end

end


# α_to_Lxᵅ(Y, index_array,k,mc) =  map(γ -> Y[[k,γ[mc[k]]]], index_array)
# """A ∈ (ℕⁿ)ᵃˣᵇ, B ∈ (ℕⁿ)ᶜˣᵈ --> D ∈ (ℕⁿ)ᵃᶜˣᵇᵈ : D₍ᵢⱼ,ₖₕ₎ = Aᵢₖ + Bⱼₕ"""
# function expo_kron(A,B)
#     n₁,n₂ = size(A)
#     D = [B + repeat( [A[i,j]] , inner = (1,1), outer = size(B)) for i in 1:n₁ , j in 1:n₂ ]
#     return cat([cat(D[i,:]...,dims=2) for i in 1:n₁]...,dims=1)
# end

# ################ Hannkel constraints for sparsity
# """For each monomial xᶜ ∈ [x]ₜ[x]ₜᵀ that occures more than twice (accounting for symmetry)
#    We generate the set of affine expressions: {∑ₖ Yᵏₐᵦ | a + b = c}"""
# function make_Hankel_con(A,t,Y) 
#     Yᵏs = moments.make_Yᵏs(t,A)
#     Hank_Yks = get_Hank_Yᵏs(t,A)
#     return[[sum(Y[pena(e,Yᵏs)]) for e in eq] for eq in Hank_Yks]
# end
# pena(ent,Yᵏs) = map(αβ -> Yᵏs[αβ[1]][αβ[2]],ent)
# """Returns a matrix of the size of the moment matrix L([x]ₜ[x]ₜᵀ)[I] but with entries
# integers. Entries with the same integers have the same monomial"""
# function get_Hank_grid(t,M)
#     meₜₜ = moments.make_mon_expo((t,t),M)
#     mm_supp = moments.get_mom_mat_supp(t,M)

#     γ_set = unique(meₜₜ)
#     γ_dict = Dict(zip(γ_set,1:length(γ_set)))

#     n_m = size(meₜₜ)[1]
#     up_tri = [i ≤ j ? 1 : 0 for i ∈ 1:n_m, j ∈ 1:n_m ]
#     return [ γ_dict[γ] for γ ∈ meₜₜ] .* mm_supp .* up_tri 
# end
# """ select only monomials that occure more than once and are supported"""
# get_Hank_list(Hank_grid) = [findall(Hank_grid .== a) for a ∈ sort(unique(Hank_grid))[2:end] if length(findall(Hank_grid .== a)) > 1 ] 
# """ find the coordinates of monomials as they occure in Yᵏ"""
# function get_Hank_Yᵏs(t,M)
#     meₜ = moments.make_mon_expo(t,M)
#     Y = [sort([α,β]) for α ∈ meₜ, β ∈ meₜ]
#     Yᵏs = moments.make_Yᵏs(t,M)
#     Hank_grid = get_Hank_grid(t,M)
#     Hank_list = get_Hank_list(Hank_grid)
#     return [[map(α -> α[1],find_Yᵏ_coords(Y[co],Yᵏs)) for co ∈ h] for h ∈ Hank_list]
# end
# function find_Yᵏ_coords(αβ,Yᵏs)
#     p = length(Yᵏs)
#     find_k_s(αβ,Yᵏs,k) = findall([αβ] .==  get_Yᵏ_ind(Yᵏs[k]))
#     return [ [[k,f] for f in find_k_s(αβ,Yᵏs,k)] for k ∈ 1:p if find_k_s(αβ,Yᵏs,k) != []]  
# end
# function get_Yᵏ_ind(Yᵏ)
#     n_Y = size(Yᵏ)[1]
#     return [ i ≥ j ? Yᵏ[i,j][2] : [] for j ∈ 1:n_Y, i ∈ 1:n_Y]
# end
# #################33