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



""" Simulation """
function simulate!(x0,p,tspan)
prob = ODEProblem(model!,x0,tspan,p)
sol = solve(prob,Tsit5(),callback=model_cb)#CallbackSet(model_cb,fsw))
return convertSol(sol.t,sol.u)
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
