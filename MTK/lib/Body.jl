include("Core.jl")

struct BodyStates <: States   
    θ::Vector{Float64} #Euler angles
    ω::Vector{Float64} #Body rate    
    #H::Vector{Float64} #Angular momentum
    #T::Vector{Float64} #Torque
end

struct BodyParameters <: Parameters
    J::Matrix{Float64} #Inertia tensor for the body
    invJ::Matrix{Float64} #inverse of J
    BodyParameters(J) = new(J,inv(J))
end

struct Body <: Bodies
    x :: BodyStates 
    p :: BodyParameters
    f :: Function    
end

function BodyDynamics(dx,x,t,p)
    ω = @view x[:ω]
    J = @view p[:J]
    invJ = @view p[:invJ]
    #H = J*ω# + Hi
    dx.θ = ω
    #dx.ω = invJ*(T - ω×(J*ω)) 
    dx.ω = invJ*( - ω×(J*ω)) 
    return dx
end