include("Core.jl")

struct ReactionWheelStates <: State
    u::Float64 #Current Command
    ω::Float64 #Wheel speed about rotation axis
    Tw::Float64 #Torque of the wheel in wheel frame
    H::Vector{Float64} #Angular momentum in parent frame
    T::Vector{Float64} #Reaction torque in parent frame    
end

struct ReactionWheelParameters <: Parameter
    kt::Float64 #Current to Torque conversion
    J::Float64 #Wheel moment of inertia about rotation axis
    a::Vector{Float64} #axis of rotation in parent frame
end

struct ReactionWheel <: Actuator
    x :: ReactionWheelStates
    p :: ReactionWheelParameters
    f :: Function    
    save
end

function ReactionWheelDynamics(dx,x,p)
    u,ω = @view x[[:u,:ω]]
    kt,J,a = @view p[[:kt,:J,:a]]
    
    dxω = Tm/J    
    T = Tm*a
    H = J*ωs*a    
    return dxω,T,H
end
#=
saved_values = SavedValues(Float64,ReactionWheelStates)
save_func = (u,t,integrator) -> copy(integrator.x) # need to copy or all saved vals will be the end vals

save = SavingCallback(save_func,saved_values)
=#