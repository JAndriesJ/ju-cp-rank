using  Plots
## ------- having fun wiht chaos

function logistic_fun(x_0::Float64,r,n) 
    w = [x_0]
    for i ∈ 1:n
      x_next = r*w[end]*(1-w[end])
       push!(w, x_next)
    end
    return w
end

x_0 = 0.21
r = 3.9 # 2.7 3.1 3.5 3.9 
n = 20
w = logistic_fun(x_0,r,n)
plot(w)

u = logistic_fun(x_0,2.7,n)
plot(u)

v = logistic_fun(x_0,3.1,n)
plot(v)



y_0 =[0.21,0.47]
v = logistic_fun(y_0, 3.8, 100000)

function logistic_fun(x_0,r,n) 
    w = [x_0]
    for i ∈ 1:n
      x_next = r*w[end].*(1 .- w[end])
       push!(w, x_next)
    end
    return hcat(w...)'
end
scatter(v[:,1],v[:,2],alpha=0.1)