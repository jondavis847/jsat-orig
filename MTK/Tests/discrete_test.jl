using ModelingToolkit, DifferentialEquations

@variables t x(t)=0 y(t)=0
D = Differential(t)
U = DiscreteUpdate(t; dt = 1)

eqs = [D(x) ~ 1
       D(y) ~ 0
       y ~ Sample(t,1)(x)]

@named de = ODESystem(eqs,t)
sys = structural_simplify(de)
prob = ODEProblem(sys,[],(0,10))
sol = solve(prob)