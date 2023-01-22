using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations, StaticArrays,Dates
includet("utils.jl")
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


""" True Dynamics """
function eom!(dx, x, p, t)
    bodyTranslation!(dx, x, p, t)
    bodyRotation!(dx, x, p, t)
end

function eom_cb!(integrator)
    bodyRotation_cb!(integrator)
    bodyTranslation_cb!(integrator)
end

function time_cb!(S)    
    S.u.orbit.epoch = S.u.orbit.epoch + (S.t-S.tprev)/86400.0
    #S.u.orbit.time = Dates.julian2datetime(S.u.orbit.epoch) #datetime type as a state is causing issues
end

# Body Translation

function bodyTranslation!(dx, x, p, t)
    dx.body.r_eci = x.body.v_eci
    #dx.body.v = -p.environments.gravity.μ / (norm(x.body.r)^3) * x.body.r + x.environments.gravity.a
    dx.body.v_eci = x.environments.gravity.a
end

function bodyTranslation_cb!(S)
    S.u.body.eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), S.u.orbit.epoch, S.p.environments.geomagnetism.eop_IAU1980)    
    S.u.body.r_ecef = S.u.body.eci_to_ecef * S.u.body.r_eci
    S.u.body.lla = SVector{3}(ecef_to_geodetic(S.u.body.r_ecef))

   # S.u.body.v_b = qvrot(S.u.body.q,S.u.body.v_eci)
end

# Body Rotation 

function bodyRotation!(dx, x, p, t)
    q = x.body.q
    # Crassidis/Markley, Eq 3.20 wrt Eq 2.88
    Q = @SMatrix [
        q[4] -q[3] q[2]
        q[3] q[4] -q[1]
        -q[2] q[1] q[4]
        -q[1] -q[2] -q[3]
    ]

    dx.body.q = 0.5 * Q * x.body.ω
    dx.body.Hb = x.body.Te - x.body.Ti - cross(x.body.ω, x.body.Hs)
    #dx.body.ω = p.body.invJ*(x.body.Te - x.body.Ti - cross(x.body.ω,x.body.H))         
end

function bodyRotation_cb!(S)
    #S.u.body.H = S.p.body.J*S.u.body.ω + S.u.body.Hi     
    S.u.body.Hs = S.u.body.Hb + S.u.body.Hi
    S.u.body.ω = S.p.body.invJ * S.u.body.Hb
    if S.u.body.q[4] < 0
        S.u.body.q = -S.u.body.q
    end

    S.u.body.Te = S.u.actuators.mtb.Tb
    S.u.body.Ti = sum(S.u.actuators.rw.Tb,dims=2)
    S.u.body.Hi = sum(S.u.actuators.rw.Hb,dims=2)
end

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
