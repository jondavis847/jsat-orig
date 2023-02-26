function mcplot(sol, f1, ind=nothing)    
    tmp_t = map(x -> x.t, sol.dispersed)    
    tmp_d = [f1(sol.dispersed[i]) for i in 1:length(sol.dispersed)]
    
    if !isnothing(ind)    
        f2(x) = map(y -> y[ind], x)    
        tmp_d = f2.(tmp_d)
        lines = [scatter(; x=sol.nominal.t, y=f2(f1(sol.nominal)), mode="lines", line=attr(color=:red))]
    else
        lines = [scatter(; x=sol.nominal.t, y=f1(sol.nominal), mode="lines", line=attr(color=:red))]        
    end
    for i in 1:length(sol.dispersed)
        pushfirst!(lines, scatter(; x=tmp_t[i], y=tmp_d[i], mode="lines", opacity=0.3, line=attr(color=:grey)))
    end    
    plot(lines)
end

function jplot(t, data)
    if data isa Vector{Vector{Float64}}
        lines = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 0)
        for i in axes(data[1])[1]
            push!(lines, scatter(; x=t, y=map(x -> x[i], data), mode="lines"))
        end
    end
    plot(lines)
    #return lines
end

function jplot3(data, noplot=false)
    if (data isa Vector{Vector{Float64}}) || (data isa Vector{SVector{3,Float64}})
        p = scatter(; x=map(x -> x[1], data), y=map(x -> x[2], data), z=map(x -> x[3], data), type="scatter3d", mode="lines")
        plot(p)
    end
end
