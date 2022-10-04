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
        Ti - internal torques
        Hi - internal momentums
        TODO: make mass properties separate and an input, needed for İω term
    outputs:
"""
function RigidBody(;i=1000.0*I(3), θ0 = zeros(3), ω0 = zeros(3), Ti0 = zeros(3), Hi0 = zeros(3), name)
    @named input_T = RealInput(u_start = zeros(3), nin = 3)
    @named input_Ti = RealInput(u_start = Ti0, nin = 3)
    @named input_Hi = RealInput(u_start = Hi0,nin = 3)

    i = scalarize(only(@parameters i[1:3,1:3] = i))
    θ, ω = scalarize.(@variables(θ(t)[1:3] = θ0, ω(t)[1:3] = ω0))

    T = scalarize(input_T.u)
    Ti = scalarize(input_Ti.u)
    Hi = scalarize(input_Hi.u)

    H = i*ω + Hi
    eqs = [
        D.(θ) .~ ω
        D.(ω) .~ inv(i) * (T - Ti - ω×H)
    ]
    return compose(ODESystem(eqs; name = name), input_T,input_Ti,input_Hi)
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
function ReactionWheel(;J = 0.25 ,kt = 1,θ = I(3),ωs0 = 100*2*pi/60,name)
#0.231 kg m^2
    unit_vec = θ*[1,0,0]
    @named input_u = RealInput(u_start = 0, nin = 1) 
    @named output_T = RealOutput(u_start = zeros(3), nout = 3)
    @named output_H = RealOutput(u_start = zeros(3), nout = 3)
    @variables ωs(t) = ωs0
    kt,J,θ = scalarize.(@parameters kt=kt J=J θ=θ)

    current_cmd = input_u.u
    Tm = kt*current_cmd
    T = scalarize(Tm*unit_vec)
    H = scalarize(J*ωs*unit_vec)
    eqs = [
        D.(ωs) ~ Tm/J    
        scalarize(output_T.u) .~ T
        scalarize(output_H.u) .~ H
    ]
    return compose(ODESystem(eqs, name=name), output_T, output_H, input_u)
end

mutable struct RW    
    J::Float64
    kt::Float64
    a::Vector{Float64}
    ωs0::Float64
end

function make(component::Vector{RW})
    n = length(component)
    @named input_u = RealInput(u_start = zeros(n), nin = n) 
    @named output_T = RealOutput(u_start = zeros(3), nout = 3)
    @named output_H = RealOutput(u_start = zeros(3), nout = 3)
    @variables ωs(t)[1:n] = getfield.(component,:ωs0)
    kt,J = scalarize.(@parameters kt[1:n]=getfield.(component,:kt) J[1:n]=getfield.(component,:J))
    
    a_val = reduce(hcat,getfield.(component,:a))
    for i in 1:n a_val[:,i] = normalize(a_val[:,i]) end
    a = scalarize(only(@parameters a[1:3,1:n] = a_val))   

    T_tmp = []
    H_tmp = []
    Tm = []
    #current_cmd = input_u.u
    for i in 1:n
        #print(kt)
        #print(current_cmd)
        Tm = kt[i]*input_u.u[i]
        push!(T_tmp,Tm*a[:,i])
        push!(H_tmp,J[i]*ωs[i]*a[:,i])
    end
    print(size(T_tmp))
    if n == 1
        T = scalarize(T_tmp)
        H = scalarize(H_tmp)
    else
        T = scalarize(sum(T_tmp))
        H = scalarize(sum(H_tmp))
    end

    eqs = [
        D.(ωs) .~ Tm./J    
        scalarize(output_T.u) .~ T
        scalarize(output_H.u) .~ H
    ]
    return compose(ODESystem(eqs, name=:ReactionWheels), output_T, output_H, input_u)
end

function make(component::RW)
    make([component])
end

function RwCommand(times,durations,values;name)
    n = 3
    @named output_u = RealOutput(u_start = zeros(n), nout = n)
    times,durations,values = @parameters times[1:n] = times durations[1:n] = durations values[1:n] = values
    eqs = []
    for i in 1:n        
        push!(eqs,ifelse((t > times[i]) & (t < (times[i] + durations[i])), [output_u.u[i] ~ values[i]], [output_u.u[i] ~ 0]))                
    end
    return compose(ODESystem(eqs, name=name), output_u)
end
    """
Test function to verify results
"""

function test(sys)
    prob = ODEProblem(sys, Pair[], (0.0,10.0))
    sol = solve(prob)
    return sol
end

"""
Make connections to make sys
"""
@named rb = RigidBody()

@named thr = Thruster()
@named thr_command = Step(start_time = 20.0, duration = 1.0)
@named thr_sys = ODESystem([
        connect(thr_command.output,thr.input_u)
        ],
    t, systems = [thr, thr_command])
#thr_sys_s = structural_simplify(thr_sys)

@named rw_command = Step(start_time = 3.0, duration = 1.0)
@named rw = ReactionWheel()
@named rw_sys = ODESystem([
        connect(rw_command.output,rw.input_u)    
        ],
    t, systems = [rw, rw_command])
#rw_sys_s = structural_simplify(rw_sys)

@named sys = ODESystem([
    connect(rw_sys.rw.output_H,rb.input_Hi)
    connect(rw_sys.rw.output_T,rb.input_Ti)
    connect(thr_sys.thr.output_T,rb.input_T)
], t, systems = [rw_sys, thr_sys, rb])

sys_s = structural_simplify(sys)

"""
TBD
"""
#inertia = [
#    669.42  65.120  15.613
#    65.120  1431.23 -39.827
#    15.613  -39.827 1551.76 
#]