using ComponentArrays, DifferentialEquations, LinearAlgebra

J = 1000.0*I(3)
kp = 1000
kd = 1000
θr = zeros(3)
ωr = zeros(3)
function eom!(dx,x,p,t)
    θ = @view(x[:θ])
    ω = @view(x[:ω])        
    dx = rotation(dx,x,p.T)             
end
function fsw!(integrator)
    integrator.p.T = controller(integrator.u)
end
function controller(x)
    θ = @view(x[:θ])
    ω = @view(x[:ω])  
    T = kp*(θr-θ) + kd*(ωr-ω)    
    return T
end
function rotation(dx,x,T)
    θ = @view(x[:θ])
    ω = @view(x[:ω])
    
    H = J*ω #+ Hi
    dx.θ = ω
    dx.ω = inv(J)*(T - ω×H) 
    return dx
end
x0 = ComponentArray(
    θ = zeros(3),
    ω = [pi/4,0,0],    
    )

p = ComponentArray(T=zeros(3))
prob = ODEProblem(eom!,x0,(0,80),p,)

cb = PeriodicCallback(fsw!,1)
sol = solve(prob,Tsit5(),callback = cb)