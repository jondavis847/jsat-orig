using ModelingToolkit,DifferentialEquations,IfElse

@variables t
x = @variables u(t)=0. [irreducible=true]

eq = [ u ~ IfElse.ifelse(t > 1., 1., 0.)]
@named s = ODESystem(eq,t,x,[])
_s = structural_simplify(s)

prob = ODEProblem(s,[],(0.,10.))
sol = solve(prob)