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
    @named input_T = RealInput(u_start = zeros(3), nin = 3)
    i = scalarize(only(@parameters i[1:3,1:3] = i))
    θ, ω = scalarize.(@variables(θ(t)[1:3] = θ0, ω(t)[1:3] = ω0))

    T = scalarize(input_T.u)
    eqs = [
        D.(θ) .~ ω
        D.(ω) .~ inv(i) * (T - ω×(i*ω))
    ]
    return compose(ODESystem(eqs; name = name), input_T)
end

"""
Thruster 
    states:
    parameters:
        F - magnitude of the thruster force
        R - position vector to thruster location in reference frame
        θ - DCM from thruster line of action to reference frame (TODO: maybe use rotations.jl one day)
    inputs:
        u - thruster command (0 = off, 1 = on)
    outputs:
        T - resultant external torque in reference frame
"""
function Thruster(F=100.0, R=[0.0,5.0,0.0], θ=1.0*I(3);name) 
    @named output_T = RealOutput(u_start = zeros(3), nout = 3)  
    @named input_u = RealInput() 
    F, R, θ = scalarize.(@parameters F=F R[1:3]=R θ[1:3,1:3] = θ)

    F̄ = [F, 0, 0]
    eqs = scalarize(output_T.u) .~ scalarize(input_u.u * (R × (θ*F̄)))
    return compose(ODESystem(eqs, name=name), output_T, input_u)
end

"""
ReactionWheel   
    states:
        ωs - wheel speed

    parameters:
        kt - motor torque constant
        r - location to wheel cg
        θ - rotation from wheel frame to reference frame, wheel assumed to spin about wheel frame x axis
        i - wheel inertia
        TODO: Wheel equations, what parameters do I need?
    inputs:
        input_u - current
    outputs:
        h - internal momentum in the reference frame
        ḣ - internal torque
"""
function ReactionWheel(;i = 0.25 ,kt = 1,θ = I(3),ωs0 = 100*2*pi/60,name)
#0.231 kg m^2
    unit_vec = θ*[1,0,0]
    h0 = i*ωs0.*unit_vec
    @named input_u = RealInput() 
    @named output_T = RealOutput(u_start = zeros(3), nout = 3)
    @named output_h = RealOutput(u_start = zeros(3), nout = 3)
    ωs, h = scalarize.(@variables(ωs(t) = ωs0, h(t)[1:3] = h0))
    kt,i,θ = scalarize.(@parameters kt=kt i=i θ=θ)

    Tm = kt*input_u.u
    T = Tm*unit_vec
    h = i*ωs*unit_vec
    eqs = [
        D.(ωs) ~ Tm/i    
        scalarize(output_T.u) .~ scalarize(T)
        scalarize(output_h.u) .~ scalarize(h)
    ]
    return compose(ODESystem(eqs, name=name), output_T, output_h, input_u)
end
"""
Make connections to make sys
"""
@named rb = RigidBody()
@named thr = Thruster()
@named thr_command = Step(start_time = 3.0, duration = 1.0)
@named thr_sys = ODESystem([
        connect(rb.input_T, thr.output_T)
        connect(thr_command.output,thr.input_u)],
    t,
    systems = [rb, thr, thr_command],
)
thr_sys_s = structural_simplify(thr_sys)

@named rw_command = Step(start_time = 3.0, duration = 1.0)
@named rw = ReactionWheel()
@named rw_sys = ODESystem([
        connect(rw_command.output,rw.input_u)],
    t,
    systems = [rw, rw_command],
)

rw_sys_s = structural_simplify(rw_sys)


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
#inertia = [
#    669.42  65.120  15.613
#    65.120  1431.23 -39.827
#    15.613  -39.827 1551.76 
#]