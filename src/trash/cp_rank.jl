module cp_rank

proj_dir = dirname(@__FILE__)*"\\"
include(proj_dir*"cp_model.jl")

# using .cp_model  

export get_ξₜᶜᵖ   
## Run computations
function get_ξₜᶜᵖ(M,t,flavour) 
    contains(flavour,"dag") ? dag=true : dag=false
    contains(flavour,"xx") ? xx=true : xx=false
    contains(flavour,"G") ? G_con=true : G_con=false
    if contains(flavour,"id")
        mod = cp_model.modelξₜᶜᵖⁱᵈ(M,t,G_con=G_con,dag=dag,xx=xx)
    elseif  contains(flavour,"wsp")
        mod = cp_model.modelξₜᶜᵖˢᵖ(M,t,G_con=G_con,dag=dag,xx=xx,isWeak=true)
    elseif  contains(flavour,"sp")
        mod = cp_model.modelξₜᶜᵖˢᵖ(M,t,G_con=G_con,dag=dag,xx=xx)
    else
        error("Incorrect model specification")
    end
    return cp_model.computeξₜᶜᵖ(mod)
end  
function get_ξₜᶜᵖ(M ,t, flavour, save_path::String)
    ξₜᶜᵖ, s = capture_solver_output(get_ξₜᶜᵖ, (M ,t, flavour))
    write_solver_output(s,save_path)
    ex_moments = cp_model.extract_moments(ξₜᶜᵖ)
    return ξₜᶜᵖ, ex_moments
end

function capture_solver_output(func,args)
    original_stdout = stdout;
    (rd, _) = redirect_stdout();
    ξₜᶜᵖ = func(args...)
    s = []
    for rl in eachline(rd)
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
