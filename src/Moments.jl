module Moments

export  make_mon_expo, # Remove this one, such that only make_mon_expo_mat is used
        make_mon_expo_mat,
        make_mom_expo_keys


"""
input: n (integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of integers)
coments: Has a quadratic loss in runtime because for degree t we compute degree 1 t times...
"""
function make_mon_expo_arr(n::Int,t::Int, isLeq::Bool = true)
    if t == 0 # [x]₌₀
        return zeros(Int32,1,n)
    else # [x]₌ₜ
        temp = make_mon_expo_arr(n,t-1,isLeq)
        e₁ = hcat(1, zeros(Int32,1,n-1))
        output = e₁ .+ temp
        for i = 1:n-1
            output = vcat(output,circshift(e₁,(0,i)) .+ temp)
        end

        if isLeq # [x]≦ₜ
            output = vcat(temp, output)
        end
        return unique(output,dims=1)
    end
end

"""
input: n (integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ of [x]₌ₜ (array of arrays of integers)
"""
function make_mon_expo(n::Int,t::Int, isLeq::Bool = true)
    mon_expo_arr = make_mon_expo_arr(n,t,isLeq)
    mon_expo     = [r  for r in  eachrow(mon_expo_arr)]
    return mon_expo
end


"""
input: n(integer),t(integer)
output: exponents α ∈ Nⁿₜ of [x]≦ₜ[x]≦ₜᵀ or [x]₌ₜ[x]₌ₜᵀ where , x = (x₁,x₂,...,xₙ)
"""
function make_mon_expo_mat(n::Int,t::Tuple{Int64,Int64},isLeq::Bool = true)
    mon1      = make_mon_expo(n,t[1], isLeq)
    mon2      = make_mon_expo(n,t[2], isLeq)
    nb_mon1   = length(mon1)
    nb_mon2   = length(mon2)
    xxᵀₜ_vec = [ mon1[i]  + mon2[j] for i in 1:nb_mon1 for j in 1:nb_mon2]

    xxᵀₜ     = reshape(xxᵀₜ_vec, (nb_mon1, nb_mon2) )
    return xxᵀₜ
end

"""
input:  n(integer),t(integer)
output: array: unique exponents in [x]≦ₜ[x]≦ₜᵀ γ ∈ N_2t^n, values are indeces in
                    moment matrix array of (α,β) ∈ (N_2t^n)^2 such that α + β = γ
"""
function make_mom_expo_keys(n::Int,t::Int)
    mon_vec = make_mon_expo(n,t)
    return unique( [ α + β for α in mon_vec, β in mon_vec ])
end


end
