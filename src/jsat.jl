using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations, StaticArrays,Dates
includet("utils.jl")
includet("dynamics.jl")
includet("environments.jl")
includet("orbit.jl")
includet("actuators.jl")
includet("fsw.jl")

lograte = 0.1

""" Model """
function model!(dx, x, p, t)
    #environments!(dx,x,p,t)
    actuators!(dx, x, p, t)
    eom!(dx, x, p, t)
end

function model_cb!(integrator)
    time_cb!(integrator)
    eom_cb!(integrator)
    environments_cb!(integrator)
    actuators_cb!(integrator)   # make sure act is after env (mtb) 
    fsw!(integrator)
end

model_cb = PeriodicCallback(model_cb!,lograte,save_positions = (true,false))

""" Initialize """
function initModelParams(in,MC=false)        
    d = Dict() #initialize output
    v(x) = !MC ? x.value : rand(x.dist) #function for pulling value or rand from dist
    
    for k in keys(in)        
        if in[k] isa NamedTuple
            # if it's another NT, dig deeper
            merge!(d,Dict(k => initModelParams(in[k],MC)))        
        elseif in[k] isa ModelParameter
            # if it's just an MP, get the val
            merge!(d,Dict(k=> v(in[k])))                
        elseif (in[k] isa Vector{ModelParameter})
            # if it's a vector, loop 1d
            merge!(d,Dict(k => [v(in[k][j]) for j in eachindex(in[k])]))
        elseif (in[k] isa Matrix{ModelParameter})
            # if it's a matrix, loop 2d
            tmp = zeros(size(in[k]))
            for i in axes(tmp)[1]
                for j in axes(tmp)[2]
                    tmp[i,j] = v(in[k][i,j])
                end
            end
            merge!(d,Dict(k => tmp))                
        else
            #if it's something else, just put it back
            merge!(d,Dict(k => in[k]))
        end
    end    
    return NamedTuple(d) #convert back to an NT
end


""" Simulation """
function simulate!(x0,p,tspan; nruns = 1)    
    p0 = initModelParams(p,false)
    prob = ODEProblem(model!,x0,tspan,p0)
    sol = solve(prob,Tsit5(),callback=model_cb)
    sol = convertSol(sol.t,sol.u)
    
    sols = Vector{typeof(sol)}(undef,nruns)
    if nruns > 1        
        ps = Vector{typeof(p0)}(undef,nruns)
        for n in 1:nruns
            ps[n] = initModelParams(p,true)
            new_prob = remake(prob,p=ps[n])
            new_sol = solve(new_prob,Tsit5(),callback=model_cb)
            sols[n] = convertSol(new_sol.t,new_sol.u)
        end
    end
    return pushfirst!(sols,sol)
end

function recursive_sim!(x0,p,tspan,interval)
    n = (tspan[2]-tspan[1])/interval
    t = [tspan[1]]
    u = [x0]
    prob = ODEProblem(model!,x0,tspan,p)
    t0 = tspan[1]
    tf = t0+interval
    uf = x0
    for i in 1:n
        prob = remake(prob,u0 = uf, tspan = (t0,tf))
        sol = solve(prob,Tsit5(),callback=model_cb)
        t0 = sol.t[end]
        tf = t0+interval
        push!(t,t0)
        uf = sol.u[end]
        push!(u,uf)
    end
    return convertSol(t,u)
end

function convertSol(t,u)
    u_out = digdeeper(u)
    return (t=t,u=u_out)
end
function digdeeper(in)
    d = Dict()
    for k in keys(in[1])
        if in[1][k] isa ComponentVector
           merge!(d,Dict(k => digdeeper(map(x->x[k],in))))
        else
            merge!(d,Dict(k=> map(x->x[k],in)))
        end
    end
    return NamedTuple(d)
end

