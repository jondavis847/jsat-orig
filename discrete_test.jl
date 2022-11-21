using ModelingToolkit, DifferentialEquations

@parameters t a b c d
@variables x(t) y(t)
δ = Differential(t)
Δ = Difference(t; dt = 0.1)
U = DiscreteUpdate(t; dt = 0.1)
eqs = [δ(x) ~ a * x - b * x * y
       δ(y) ~ -c * y + d * x * y
       Δ(x) ~ y
       U(y) ~ x + 1]
@named de = ODESystem(eqs, t, [x, y], [a, b, c, d])
