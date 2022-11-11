using ComponentArrays, DifferentialEquations, LinearAlgebra
mutable struct RigidBody
    J::Matrix{Float64}
    θ::Vector{Float64}
    ω::Vector{Float64}
end
    
θr = zeros(3)
ωr = zeros(3)
#kp::Float64 = 1000
 #   kd::Float64 = 1000
function eom!(dx,x,p,t)
    θ = @view(x[:θ])
    ω = @view(x[:ω])        
    dx = rotation(dx,x,p.T)             
end

function rotation(dx,x,T)
    θ = @view(x[:θ])
    ω = @view(x[:ω])
    
    H = J*ω #+ Hi
    dx.θ = ω
    dx.ω = inv(J)*(T - ω×H) 
    return dx
end
