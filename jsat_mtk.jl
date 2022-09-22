using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks, DifferentialEquations, LinearAlgebra, Rotations, Plots
using Symbolics: scalarize 

@variables t
D = Differential(t) 

"""
RigidBody
    states:
        θ - euler angles (TODO: make this quat one day)
        ω - angular velocity
    parameters:
        i - inertia tensor (little i to not conflict with identity matrix function)
        θ0 - initial euler angles
        ω0 - initial angular velocity
    inputs:
        T - external torques
        h - internal momentums
        TODO: make mass properties separate and an input, needed for İω term
    outputs:
"""
function RigidBody(;i=1000.0*I(3), θ0 = zeros(3), ω0 = zeros(3), name)
    @named input_T = RealInput(u_start = zeros(3))
    @parameters i[1:3,1:3]
    @variables θ(t)[1:3] = θ0
    @variables ω(t)[1:3] = ω0

    T = input_T.u
    eqs = [
        scalarize(D.(θ) .~ ω)
        scalarize(D.(ω) .~ inv(scalarize(i))*(scalarize(T).-scalarize(cross(ω,i*ω)))) #gotta be a better way than scalarizing everything independently, but this works for now
        # D.(ω) .~ inv(I)*(T-h-İ*ω-cross(ω,I*ω))
    ]
    return compose(ODESystem(eqs, t, [θ...;ω...], [i...];name = name), input_T)
end

"""
Thruster 
    states:
    parameters:
        Fmag - magnitude of the thruster force
        R - position vector to thruster location in reference frame
        rotation - DCM from thruster line of action to reference frame (TODO: maybe use rotations.jl one day)
    inputs:
        u - thruster command (0 = off, 1 = on)
    outputs:
        T - resultant external torque in reference frame
"""
function Thruster(Fmag=1.0,R=[1.0,0,0],rotation=1.0*I(3);name) 
    @named output_T = RealOutput(u_start = zeros(3))  
    @named input_u = RealInput() 
    @parameters Fmag, R[1:3], F[1:3], rotation[1:3,1:3]

    eqs = scalarize(output_T.u .~ input_u.u*cross(R,rotation*[Fmag,0,0]))
    compose(ODESystem(eqs,t,[], [Fmag,R...,F...,rotation...],name=name),output_T,input_u)
end

"""
Make connections to make sys
"""
@named rb = RigidBody()
@named th = Thruster()
@named command = Step(start_time = 3.0, duration = 1.0)
@named sys = ODESystem([
    connect(rb.input_T, th.output_T)
    connect(command.output,th.input_u)
    ],t,systems = [rb,th])


"""
Test function to verify results
"""
function test(sys)
prob = ODEProblem(sys, Pair[], (0.0,10.0))
sol = solve(prob)
return sol
end

"""
TBD
"""
#thrusters_on = [3.0] => [T[1] ~ 100]
#thrusters_off = [4.0] => [T[1] ~ 0]
#thrusters = [thrusters_on, thrusters_off]




#inertia = [
#    669.42  65.120  15.613
#    65.120  1431.23 -39.827
#    15.613  -39.827 1551.76 
#]