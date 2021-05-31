module Constraints

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"Utility.jl")
include(proj_dir*"Moments.jl")
using .Utility
using .Moments

export  make_dag_con,
        make_loc_con,
        make_xx_con,
        make_weakG_con,
        make_G_con


"""
input: A(Data matrix),t(integer),Lx(JuMP variable)
output: dictionary: keys:
comment: L(gu) ≥ 0 for g ∈ {1} ∪ Sᶜᵖ_A and u ∈ [x]_2t−deg(g)
"""
function make_dag_con(A,t,Lx)
    n = size(A)[1]
    dag_con = Dict()

    mom₂ₜ     = Moments.make_mon_expo(n,2*t)
    nb2t      = lastindex(mom₂ₜ)
    dag_con[(0,0)] = [ 1*Lx[mom₂ₜ[i] ] for i in  1:nb2t] # This is the non negativity of the moments i.e. products of x's

    deg_g = 2
    mom₂ₜ₋₂    = Moments.make_mon_expo(n,2*t-deg_g)
    nb2t₋₂ = lastindex(mom₂ₜ₋₂)
    @assert nb2t₋₂ > 1 "Contraints do not exist for t =1!"
    for k in 1:n
        eₖ = Utility.eₖ(n,k)
        # Dagger constraints: L((√Aₖₖ xₖ - xₖ²)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
        sqrMₖₖ = sqrt(A[k,k])
        dag_con[(k,k)] = [ sqrMₖₖ*Lx[eₖ + mom₂ₜ₋₂[i]] - Lx[2*eₖ + mom₂ₜ₋₂[i]] for i in  1:nb2t₋₂ ]
        for h in (k+1):n
            eₕ = Utility.eₖ(n,h)
            # Dagger constraints: L((Aₖₕ  - xₖxₕ)⋅u) ≧ 0 for u ∈ [x]₂ₜ₋₂
            dag_con[(k,h)] = [ A[k,h]*Lx[mom₂ₜ₋₂[i]] - Lx[eₖ + eₕ + mom₂ₜ₋₂[i]]     for i in  1:nb2t₋₂]
        end
    end
    return dag_con
end

"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (i,j) ∈ [n]²
                    vals:
comment: L ≥ 0 on M₂ₜ(S^cp_A )
(M_2t-2(gL) )_αβ =   √(Aᵢᵢ) x^(γ + eᵢ)  -  x^(γ + 2*eᵢ) M_2
(M_2t-2(gL) )_αβ =   (Aᵢⱼ) x^γ   -  x^(γ + e₁ + eⱼ) """
function make_loc_con(A,t,Lx)
    n       = size(A)[1]
    LMB     = Moments.make_mon_expo(n, t - 1)
    nb_mon  = size(LMB)[1]
    loc_con = Dict()
    for k in 1:n
        eₖ = Utility.eₖ(n,k)
        # Constraint: diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((√Aₖₖ xₖ - xₖ²)⋅L)
        sqrAₖₖ = sqrt(A[k,k])
        loc_con[(k,k)] = [sqrAₖₖ*Lx[LMB[i] + LMB[j] + eₖ] - Lx[LMB[i] + LMB[j] + 2*eₖ]   for i in 1:nb_mon, j in 1:nb_mon ]
        for h in (k+1):n
            eₕ = Utility.eₖ(n,h)
            # Constraint: off diagonal of L ≧ 0 on 𝑀(Sᶜᵖ)   : M((Aₖₕ - xₖxₕ)⋅L)
            loc_con[(k,h)] = [ A[k,h]*Lx[LMB[i] + LMB[j]] - Lx[LMB[i] + LMB[j] + eₖ + eₕ] for i in 1:nb_mon,  j in 1:nb_mon ]
        end
    end
    return loc_con
end

"""
input: A(data array), LMB(moment exponent vector), Lx(JuMP variable)
output: dictionary: Keys: (h,k) ∈ [n]², h ≠ k
                    values:
"""
function make_xx_con(A,t,Lx)
    n       = size(A)[1]
    LMB     = Moments.make_mon_expo(n,t-1)
    nb_mon  = size(LMB)[1]
    xx_con  = Dict()
    for k in 1:n
        eₖ = Utility.eₖ(n,k)
        for h in (k+1):n
            eₕ = Utility.eₖ(n,h)
            # Localizing xx constraint: M(xₖxₕ⋅L)
            xx_con[(k,h)] = [ Lx[LMB[i] + LMB[j] + eₖ + eₕ] for i in 1:nb_mon, j in 1:nb_mon ]
        end
    end
    return xx_con
end

## Tensor constraints
"""
Input: A(data matrix),t(Integer),Lx(JuMP variable)
Output: L((M-([x]₌₁[x]₌₁ᵀ))⊗([x]₌ₗ[x]₌ₗᵀ)))for l ∈ 0,1,...,t-1.
= M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))
"""
function make_weakG_con(A,t,Lx)
    n = size(A)[1]
    weakG_con = Dict()
    LMBexp_1 =  make_mon_expo_mat(n,(1,1),false)
    for ℓ in 1:(t-1)
        LMBexp_ℓ          = Moments.make_mon_expo_mat(n,(ℓ,ℓ),false)   #exponents of [x]₌ₗ[x]₌ₗᵀ
        LMBexp_1ℓ         = Utility.var_kron(LMBexp_1,LMBexp_ℓ)    #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ)
        LMB_ℓ             = Utility.index_to_var(Lx,LMBexp_ℓ)      # L([x]₌ₗ[x]₌ₗᵀ)
        LMB_1ℓ            = Utility.index_to_var(Lx,LMBexp_1ℓ)     # L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ))

        A_tens_LMB_ℓ      = kron(A,LMB_ℓ)                  # M⊗L([x]₌ₗ[x]₌ₗᵀ)
        weakG_con[ℓ]      = A_tens_LMB_ℓ - LMB_1ℓ          # M⊗L([x]₌ₗ[x]₌ₗᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]₌ₗ[x]₌ₗᵀ) ⪰ 0, ℓ ∈ 0,1,t-deg(g)/2
    end
    return weakG_con
end

"""M(G ⊗ L) ⪰ 0 constraints
input: A(data matrix),t(Integer),Lx(JuMP variable)
Assumption: G = A-[x]₌₁[x]₌₁ᵀ
output: A⊗L([x]ₜ[x]ₜᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
"""
function make_G_con(A,t,Lx)
    n = size(A)[1]
    LMBexp_1           = Moments.make_mon_expo_mat(n,(1,1),false) #exponents of [x]₌₁[x]₌₁ᵀ
    LMBexpₜ₋₁          = Moments.make_mon_expo_mat(n,(t-1,t-1),true)#exponents of [x]ₜ₋₁[x]ₜ₋₁ᵀ
    LMBₜ₋₁             = Utility.index_to_var(Lx,LMBexpₜ₋₁)    #L([x]ₜ₋₁[x]ₜ₋₁ᵀ)

    LMBexp_1ₜ₋₁        = Utility.var_kron(LMBexp_1,LMBexpₜ₋₁)  #exponents of([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ)
    LMB_1ₜ₋₁           = Utility.index_to_var(Lx,LMBexp_1ₜ₋₁)   # L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))

    G_con = kron(A,LMBₜ₋₁) - LMB_1ₜ₋₁             # A⊗L([x]ₜ₋₁[x]ₜ₋₁ᵀ) - L(([x]₌₁[x]₌₁ᵀ)⊗([x]ₜ₋₁[x]ₜ₋₁ᵀ))
    return G_con
end

end
