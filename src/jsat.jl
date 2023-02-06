using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations, StaticArrays, Dates, PlotlyJS
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

includet("utils.jl")
includet("dynamics.jl")
includet("environments.jl")
includet("orbit.jl")
includet("actuators.jl")
includet("fsw\\fsw.jl")

lograte = 0.1

""" Model """
function model!(dx, x, p, t)
    #environments!(dx,x,p,t)
    actuators!(dx, x, p, t)
    eom!(dx, x, p, t)
end

function model_cb!(integrator)    
    time_cb!(integrator)
    orbit_cb!(integrator)
    eom_cb!(integrator)
    environments_cb!(integrator)
    actuators_cb!(integrator)   # make sure act is after env (mtb) 
    fsw!(integrator)
end

model_cb = PeriodicCallback(model_cb!, lograte, save_positions=(false, false))

""" Initialize """
function initModelParams(in, MC=false)
    d = Dict() #initialize output
    v(x) = !MC ? x.value : rand(x.dist) #function for pulling value or rand from dist

    for k in keys(in)
        if in[k] isa NamedTuple
            # if it's another NT, dig deeper
            merge!(d, Dict(k => initModelParams(in[k], MC)))
        elseif in[k] isa ModelParameter
            # if it's just an MP, get the val
            merge!(d, Dict(k => v(in[k])))
        elseif (in[k] isa Vector{ModelParameter})
            # if it's a vector, loop 1d
            merge!(d, Dict(k => [v(in[k][j]) for j in eachindex(in[k])]))
        elseif (in[k] isa Matrix{ModelParameter})
            # if it's a matrix, loop 2d
            tmp = zeros(size(in[k]))
            for i in axes(tmp)[1]
                for j in axes(tmp)[2]
                    tmp[i, j] = v(in[k][i, j])
                end
            end
            merge!(d, Dict(k => tmp))
        else
            #if it's something else, just put it back
            merge!(d, Dict(k => in[k]))
        end
    end
    return NamedTuple(d) #convert back to an NT
end


""" Simulation """

function simulate(x0, p, tspan; nruns=1, dt=0.1)
    p0 = initModelParams(p, false)
    prob = ODEProblem(model!, x0, tspan, p0, callback=model_cb)

    sol = solve(prob, Tsit5(), saveat=dt)
    sol = convertSol(sol.t, sol.u)
    if nruns > 1
        eprob = EnsembleProblem(prob, prob_func=prob_func, output_func=output_func)
        sols = solve(eprob, Tsit5(), trajectories=nruns, saveat=dt)
        sol_out = (nominal=sol,dispersed=sols)
    else
        sol_out = sol
    end
    return sol_out
end

function prob_func(prob, i, repeat)
    print("Simulating run $i\n")
    remake(prob, p=initModelParams(p, true))
end
function output_func(sol, i)
    sol = convertSol(sol.t, sol.u)
    return (sol, false)
end
#=
function recursive_sim(prob, tspan, interval)
    n = (tspan[2] - tspan[1]) / interval
    t = Float64[tspan[1]]
    u = [prob.u0]
    t0 = tspan[1]
    tf = t0 + interval
    uf = prob.u0
    for i in 1:n
        prob = remake(prob, u0=uf, tspan=(t0, tf))
        sol = solve(prob, Tsit5(), callback=model_cb)
        t0 = sol.t[end]
        tf = t0 + interval
        push!(t, t0)
        uf = sol.u[end]
        push!(u, uf)
    end
    return convertSol(t, u)
end
=#

function convertSol(t, u)
    u_out = digdeeper(u)
    return (t=t, u=u_out)
end
function digdeeper(in)
    d = Dict()
    for k in keys(in[1])
        if in[1][k] isa ComponentVector
            merge!(d, Dict(k => digdeeper(map(x -> x[k], in))))
        else
            merge!(d, Dict(k => map(x -> x[k], in)))
        end
    end
    return NamedTuple(d)
end


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
