includet("lib\\Body.jl")

using DifferentialEquations, LinearAlgebra
x_b = BodyStates([pi/2,0,0],[0,pi/10,0])
p_b = BodyParameters(1000*I(3))

b = Body(x_b,p_b,BodyDynamics)
sc = Spacecraft([],[b],[],[],[])

prob = ODEProblem(sc.bodies[1].f,   )
