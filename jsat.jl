using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations, StaticArrays
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
    environments_cb!(integrator)
    actuators_cb!(integrator)
    eom_cb!(integrator)
end

model_cb = PeriodicCallback(model_cb!,lograte)


""" True Dynamics """
function eom!(dx, x, p, t)
    bodyTranslation!(dx, x, p, t)
    bodyRotation!(dx, x, p, t)
end

function eom_cb!(integrator)
    bodyRotation_cb!(integrator)
end


# Body Translation

function bodyTranslation!(dx, x, p, t)
    dx.body.r = x.body.v
    dx.body.v = -p.environments.gravity.μ / (norm(x.body.r)^3) * x.body.r + x.environments.gravity.a
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
    dx.body.H = x.body.Te - x.body.Ti - cross(x.body.ω, x.body.H)
    #dx.body.ω = p.body.invJ*(x.body.Te - x.body.Ti - cross(x.body.ω,x.body.H))         
end

function bodyRotation_cb!(S)
    #S.u.body.H = S.p.body.J*S.u.body.ω + S.u.body.Hi     
    S.u.body.ω = S.p.body.invJ * (S.u.body.H - S.u.body.Hi)
    if S.u.body.q[4] < 0
        S.u.body.q = -S.u.body.q
    end
end

""" Simulation """
function simulate!(x0,p,tspan)
prob = ODEProblem(model!,x0,tspan,p)
sol = solve(prob,Tsit5(),callback=CallbackSet(model_cb,fsw))
return sol
end

