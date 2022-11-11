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

function rotation(dx,x,T)
    θ = @view(x[:θ])
    ω = @view(x[:ω])
    
    H = J*ω #+ Hi
    dx.θ = ω
    dx.ω = inv(J)*(T - ω×H) 
    return dx
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


x0 = ComponentArray(
    θ = zeros(3),
    ω = [pi/4,0,0],    
    )

p = ComponentArray(T=zeros(3))
prob = ODEProblem(eom!,x0,(0,80),p,)

saved_values = SavedValues(Float64,typeof(p))
save_func = (u,t,integrator) -> copy(integrator.p) # need to copy or all saved vals will be the end vals

save = SavingCallback(save_func,saved_values)

fsw = PeriodicCallback(fsw!,.1)
sol = solve(prob,Tsit5(),callback = CallbackSet(fsw,save))