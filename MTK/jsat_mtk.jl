using ModelingToolkit, DifferentialEquations, LinearAlgebra, Rotations, Plots, IfElse
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
    J::Matrix{Float64}
    θ0::Vector{Float64}
    ω0::Vector{Float64}    
end

function make(component::RigidBody)   
    @variables Te(t)[1:3] = zeros(3) [input = true]
    @variables Ti(t)[1:3] = zeros(3) [input = true]
    @variables Hi(t)[1:3] = zeros(3) [input = true]
    @variables θ(t)[1:3] = component.θ0 [output = true]
    @variables ω(t)[1:3] = component.ω0 [output = true]
    Te,Ti,Hi,θ,ω = scalarize.([Te,Ti,Hi,θ,ω])            

    J = scalarize(only(@parameters J[1:3,1:3] = component.I))

    H = J*ω + Hi
    eqs = [
        D.(θ) .~ ω
        D.(ω) .~ inv(J) * (Te - Ti - ω×H)        
    ]    
    return ODESystem(eqs,t,name=:rb)
end

function makeSimple(component::RigidBody)
    @variables T(t)[1:3] = zeros(3) [input = true]
    @variables θ(t)[1:3] = component.θ0
    @variables ω(t)[1:3] = component.ω0 [output = true]
    T,θ,ω = scalarize.([T,θ,ω])            
    J = scalarize(only(@parameters J[1:3,1:3] = component.J))

    H = J*ω
    eqs = [
        D.(θ) .~ ω
        D.(ω) .~ inv(J) * (T-ω×H)        
    ]    
    return ODESystem(eqs,t,name=:rb)
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
mutable struct Thruster
    F::Float64
    R::Vector{Float64}
    θ::Matrix{Float64}
end
function make(component::Thruster) 
    @variables u(t) = 0 [input = true]
    @variables T(t)[1:3] = zeros(3) [output = true]
    T = scalarize(T)

    @parameters F = component.F
    @parameters R[1:3] = component.R
    @parameters θ[1:3,1:3] = component.θ
    F,R,θ = scalarize.([F,R,θ])
    
    F̄ = [F, 0, 0]
    eqs = zeros(3) .~ T - u*(R × (θ*F̄))

    return ODESystem(eqs,t,name=:thr)
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
    @variables (u(t))[1:n] = zeros(n) [input = true] #input current command
    @variables (T(t))[1:3] = zeros(3) [output = true] #output internal reaction torque 
    @variables (H(t))[1:3] = zeros(3) [output = true] #output internal momentum
    @variables ωs(t)[1:n] = getfield.(component,:ωs0) #array of wheel speeds
    
    @parameters kt[1:n]=getfield.(component,:kt) #current to torque motor gain
    @parameters J[1:n]=getfield.(component,:J) #wheel inertia    
    
    a_val = reduce(hcat,getfield.(component,:a))
    for i in 1:n a_val[:,i] = normalize(a_val[:,i]) end
    a = scalarize(only(@parameters a[1:3,1:n] = a_val))   

    Tm = [kt[i]*u[i] for i in 1:n]    

    eqs = scalarize([
        D.(ωs) .~ Tm./J    
        T .~ sum([Tm[i]*a[:,i] for i in 1:n])
        H .~ sum([J[i]*ωs[i]*a[:,i] for i in 1:n])
    ])
    return ODESystem(eqs,t,name=:rw)
end

function make(component::ReactionWheel)    
    return make([component])
end


"""
Step function for unit testing
    n - number of step functions out put as an array (single output)
"""
function step(;name,times,durations,values)
    @assert (length(times) == length(durations)) & (length(durations) == length(values)) "times,values,durations must all be the same length"
    n = length(times) 
    @variables (u(t))[1:n]=zeros(n) [output = true, irreducible = true]        
    @parameters begin
        times[1:n] = times
        durations[1:n] = durations
        values[1:n] = values
    end       
    
    eqs = [0 ~ u[i] - IfElse.ifelse(((times[i] < t) & ((times[i] + durations[i]) > t)), values[i], 0) for i in 1:n]    
    return ODESystem(eqs,t,name=name)
end

function controller(;name,samplerate,kd=1,kp=1)
    #x = @variables θ(t)[1:3] [input = true] ω(t)[1:3]=zeros(3) [input = true] u(t)[1:3]=zeros(3) [output = true]        
    x = @variables ω(t)[1:3]=zeros(3) [input = true] u(t)[1:3]=zeros(3) [output = true]        
    p = @parameters kd=kd kp=kp
    
    ω_ref = zeros(3)
    θ_ref = zeros(3)
    U = DiscreteUpdate(t; dt = samplerate)
    eqs = zeros(3).~ U.(u) - kd*(ω_ref-ω)# + kp*(θ_ref - θ)
    
    return ODESystem(scalarize(eqs),t,[x...;],p,name=name)
end

"""
Test function to verify results
"""
function test(sys)
    prob = ODEProblem(sys, [], (0.0,10.0))
    sol = solve(prob,Rodas4(),dtmax = 0.01)
    return sol
end


"""
Example systems
"""
#=
#Make the Rigid Body
rb = make(RigidBody(1000*I(3),zeros(3),zeros(3)))

#Make the Reaction Wheels
r1 = ReactionWheel(0.25, 1, [ 1, 0, 0], 20)
r2 = ReactionWheel(0.35, 1.1, [ 0, 1, 0], 30)
r3 = ReactionWheel(0.45, 1.2, [ 0, 0, 1], 40)
rw = make([r1,r2,r3])

    #=
#Make the Thrusters
@named thr = Thruster()
@named thr_command = Step(start_time = 9, duration = 0.5)
@named thr_sys = ODESystem([
        connect(thr_command.output,thr.input_u)
        ], t, systems = [thr, thr_command])

@named fake_thrusters = nStep(times = [Inf,Inf,Inf], durations = ones(3),values = ones(3)) 

# Connect components to Rigid Body

@named sys = ODESystem([
    connect(rw_sys.ReactionWheels.output_T, rb.input_Ti)
    connect(rw_sys.ReactionWheels.output_H, rb.input_Hi)
    connect(fake_thrusters.output_u, rb.input_T)
    ],t,systems = [rw_sys,fake_thrusters,rb])

sys_s = structural_simplify(sys)
=#
=#