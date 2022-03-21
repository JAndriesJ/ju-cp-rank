module constraints

using Test
include(dirname(@__FILE__)*"\\"*"moments.jl")
using .moments ; const mom = moments

export  α_to_Lxᵅ,
        make_dag_con,
        make_loc_con,
        make_xx_con,
        make_ideal_con,
        make_G_con

"""
input: A(Data matrix),t(integer),Lx(JuMP variable)
output: dictionary: keys:
comment: L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)
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
    momₜ₋₁  = mom.make_mon_expo(n, t-1, A)
    # if_diag(k) = [sqrt(A[k,k])*Lx[α+β+mom.eᵢ(n,k)] - Lx[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    # if_off_diag(k,h) = [A[k,h]*Lx[α+β] - Lx[α+β+mom.eᵢⱼ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
    momₜ₋₁g(j)  = mom.make_mon_expo(n, t-1, A,j)
    if_diag(k) = [sqrt(A[k,k])*Lx[α+β+mom.eᵢ(n,k)] - Lx[α+β+2*mom.eᵢ(n,k)] for α ∈ momₜ₋₁g(k), β ∈ momₜ₋₁g(k)]
    if_off_diag(k,h) = [A[k,h]*Lx[α+β] - Lx[α+β+mom.eᵢⱼ(n,k,h)] for α ∈ momₜ₋₁, β ∈ momₜ₋₁]
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
    eᵢⱼs    = map( e -> mom.eᵢⱼ(n,e[1],e[2]), mom.get_nonzero_entries(A)) 
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
    LMBexpₜ₋₁   = mom.make_mon_expo(n,(t-1,t-1),A)     #exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁      = α_to_Lxᵅ(Lx,LMBexpₜ₋₁)    #L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp₁ₜ₋₁  = expo_kron(LMBexp₁,LMBexpₜ₋₁)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB₁ₜ₋₁     = α_to_Lxᵅ(Lx,LMBexp₁ₜ₋₁)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB₁ₜ₋₁             # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    return G_con
end


"""M(xᵢxⱼL) = 0  for all {i,j} s.t. Mᵢⱼ"""
function make_ideal_con(A,t,Lx) 
    n = size(A)[1] 
    momₜ₋₁ = mom.make_mon_expo(n, t-1,A)
    eᵢⱼs = map( e -> mom.eᵢⱼ(n,e[1],e[2]), mom.get_zero_entries(A))
    return map(eᵢⱼ -> [Lx[α+β+eᵢⱼ] for α ∈ momₜ₋₁, β ∈ momₜ₋₁ ], eᵢⱼs)
end

### Utility
"""(L,[α]ᵢⱼ) → [L(xᵅ)]ᵢⱼ """
α_to_Lxᵅ(Lx, index_array) = map(α -> Lx[α],index_array)


"""A ∈ (ℕⁿ)ᵃˣᵇ, B ∈ (ℕⁿ)ᶜˣᵈ --> D ∈ (ℕⁿ)ᵃᶜˣᵇᵈ : D₍ᵢⱼ,ₖₕ₎ = Aᵢₖ + Bⱼₕ"""
function expo_kron(A,B)
    n₁,n₂ = size(A)
    D = [B + repeat( [A[i,j]] , inner = (1,1), outer = size(B)) for i in 1:n₁ , j in 1:n₂ ]
    return cat([cat(D[i,:]...,dims=2) for i in 1:n₁]...,dims=1)
end

function run_tests()
    @testset "expo_kron" begin
        n = 4
        t = (2,3)
        A = mom.make_mon_expo(n,t)
        B = mom.make_mon_expo(n,t)
        AOXB = expo_kron(A,B)
        @test size(AOXB) == size(A) .* size(B) 
    end
end

end

## Tensor constraints
# """
# Input: A(data matrix),t(Integer),Lx(JuMP variable)
# Output: L((M-([x]₌₁[x]₌₁ᵀ))⊗([x]₌ₗ[x]₌ₗᵀ)))for l ∈ 0,1,...,t-1.
# = M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))
# """
# function make_weakG_con(A,t,Lx)
#     n = size(A)[1]
#     weakG_con = Dict()
#     LMBexp_1 =  make_mon_expo(n,(1,1),isle=false)
#     for ℓ in 1:(t-1)
#         LMBexp_ℓ          = mom.make_mon_expo(n,(ℓ,ℓ),isle=false)   #exponents of [x]₌ₗ[x]₌ₗᵀ
#         LMBexp_1ℓ         = expo_kron(LMBexp_1,LMBexp_ℓ)    #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
#         LMB_ℓ             = α_to_Lxᵅ(Lx,LMBexp_ℓ)      # L([x]₌ₗ[x]₌ₗᵀ)
#         LMB_1ℓ            = α_to_Lxᵅ(Lx,LMBexp_1ℓ)     # L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))

#         A_tens_LMB_ℓ      = kron(A,LMB_ℓ)                  # M⊗L([x]₌ₗ[x]₌ₗᵀ)
#         weakG_con[ℓ]      = A_tens_LMB_ℓ - LMB_1ℓ          # M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ) ⪰ 0, ℓ ∈ 0,1,t-deg(g)/2
#     end
#     return weakG_con
# end