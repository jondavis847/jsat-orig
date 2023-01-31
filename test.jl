using DifferentialEquations

struct State
    r::Float64
    v::Float64
end

Base.length(::State) = 2
Base.oneunit(::State) = State(1,1)

x0 = State(0,1)

function f!(dx,x,p,t) 
    dx.r = x.v
    dx.v = 0.5
end

prob = ODEProblem(f!,x0,(0,10))
sol = solve(prob)