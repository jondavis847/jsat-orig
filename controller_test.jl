using ModelingToolkit,DifferentialEquations
@variables t
U = DiscreteUpdate(t,dt=1.)    
D = Differential(t)

function constant(c=1.;name)
    x = @variables u(t)=c [output = true, irreducible = true]   
    eq = D(u) ~ 0
    ODESystem(eq,t; name=name)
end

function controller(k=1.;name)   
    x = @variables y(t) [output = true,irreducible = true]  u(t) [input = true]
    p = @parameters k = k
    return ODESystem(0 ~ U(y) + k*u,t,x,p,name=name)
end

@named k = constant()
@named c = controller()

ck = compose(ODESystem(c.u ~ k.u,t,name = :ck),c,k)
#_ck = structural_simplify(ck)
#prob = ODEProblem(_ck,[],(0,10))
#sol = solve(prob)







