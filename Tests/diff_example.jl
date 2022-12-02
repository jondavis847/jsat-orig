using ModelingToolkit, DifferentialEquations
@parameters t a b c d
@variables x(t) y(t)
δ = Differential(t)
D = Difference(t; dt=0.1)
eqs = [
    δ(x) ~ a*x - b*x*y,
    δ(y) ~ -c*y + d*x*y,
    D(x) ~ y
]
@named de = ODESystem(eqs,t,[x,y],[a,b,c,d])
#@test generate_difference_cb(de) isa ModelingToolkit.DiffEqCallbacks.DiscreteCallback

# doesn't work with ODEFunction
# prob = ODEProblem(ODEFunction{false}(de),[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])

prob = ODEProblem(de,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0], check_length=false)
#@test prob.kwargs[:difference_cb] isa ModelingToolkit.DiffEqCallbacks.DiscreteCallback

sol = solve(prob, Tsit5())#; callback=prob.kwargs[:difference_cb], tstops=prob.tspan[1]:0.1:prob.tspan[2])
