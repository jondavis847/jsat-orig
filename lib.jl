using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks, DifferentialEquations, LinearAlgebra
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

mutable struct RigidBody
    I::Matrix{Float64}
    θ0::Vector{Float64}
    ω0::Vector{Float64}
    #Ti0 = zeros(3)
    #Hi0 = zeros(3)
end

function make(component::RigidBody)
    @named input_Te = RealInput(nin = 3)
    @named input_Ti = RealInput(nin = 3)
    @named input_Hi = RealInput(nin = 3)
    i = scalarize(only(@parameters i[1:3,1:3] = component.I))
    θ, ω = scalarize.(@variables(θ(t)[1:3] = component.θ0, ω(t)[1:3] = component.ω0))

    Te = scalarize(input_Te.u)
    Ti = scalarize(input_Ti.u)
    Hi = scalarize(input_Hi.u)

    H = i*ω + Hi
    eqs = [
        D.(θ) .~ ω
        D.(ω) .~ inv(i) * (Te - Ti - ω×H)
    ]
    return compose(ODESystem(eqs; name = :RigidBodys), input_Te,input_Ti,input_Hi)
end
#=function RigidBody(;i=1000.0*I(3), θ0 = zeros(3), ω0 = zeros(3), Ti0 = zeros(3), Hi0 = zeros(3), name)
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
=#
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
mutable struct Thruster
    F::Float64
    R::Vector{Float64}
    θ::Matrix{Float64}
end
function make(component::Thruster) 
    @named output_T = RealOutput(nout = 3)  
    @named input_u = RealInput(nin = 1) 
    F, R, θ = scalarize.(@parameters F=component.F R[1:3]=component.R θ[1:3,1:3] = component.θ)

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
        input_u - current command
    outputs:
        H - internal momentum in the reference frame
        T - internal torque in the reference frame
"""

mutable struct ReactionWheel    
    J::Float64
    kt::Float64
    a::Vector{Float64}
    ωs0::Float64
end

function make(component::Vector{ReactionWheel})
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
    for i in 1:n
        push!(Tm,kt[i]*input_u.u[i])
        push!(T_tmp,Tm[i]*a[:,i])
        push!(H_tmp,J[i]*ωs[i]*a[:,i])
    end
    T = sum(T_tmp)
    H = sum(H_tmp)
    eqs = scalarize([
        D.(ωs) .~ Tm./J    
        output_T.u .~ T
        output_H.u .~ H
    ])
    return compose(ODESystem(eqs, name=:rw), output_T, output_H, input_u)
end

function make(component::ReactionWheel)
    @named input_u = RealInput(u_start = 0., nin = 1) 
    @named output_T = RealOutput(u_start = zeros(3), nout = 3)
    @named output_H = RealOutput(u_start = zeros(3), nout = 3)
    @variables ωs(t) = component.ωs0
    @parameters begin
        kt = component.kt
        J = component.J
        a[1:3] = component.a
    end    
    Tm = kt*input_u.u
    T = Tm*a
    H = J*ωs*a
    
    eqs = [
        D.(ωs) ~ Tm/J    
        scalarize(output_T.u .~ T)
        scalarize(output_H.u .~ H)
    ]
    return compose(ODESystem(eqs, name=:ReactionWheels), output_T, output_H, input_u)
end

"""
Step function for unit testing
    n - number of step functions out put as an array (single output)
"""

function nStep(;name,times,durations,values)
    n = length(times)
    @named output_u = RealOutput(u_start = zeros(n), nout = n)
    @parameters begin
        times[1:n] = times
        durations[1:n] = durations
        values[1:n] = values
    end   
    eqs = []
    for i in 1:n        
        push!(eqs, 0. ~ output_u.u[i] - ifelse((t > times[i]) & (t < (times[i] + durations[i])), values[i], 0))                
    end
    return compose(ODESystem(eqs, t, name=name), [output_u])
end

#This was made to be differentiable but looks like nstep works, is simpler, and more flexible
function Step3(;name,times,durations = [1.,1.,1.])
    @named s1 = Step(start_time=times[1], duration=durations[1],smooth = true)
    @named s2 = Step(start_time=times[2], duration=durations[2],smooth = true)
    @named s3 = Step(start_time=times[3], duration=durations[3],smooth = true)
    @named output = RealOutput(nout = 3)
    eqs = scalarize(output.u .~ [s1.output.u, s2.output.u, s3.output.u])
    return compose(ODESystem(eqs, name=name),s1,s2,s3,output)  
end

"""
Test function to verify results
"""

function test(sys)
    prob = ODAEProblem(sys, Pair[], (0.0,10.0), check_length = false)
    sol = solve(prob,Rodas4())
    return sol
end
