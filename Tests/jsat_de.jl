using ComponentArrays, DifferentialEquations, LinearAlgebra, Parameters, SimulationLogs

function eom!(dx,x,p,t)
    rotation(dx,x,p,t)             
end

function rotation(dx,x,p,t)
    @unpack θ,ω = x
    J = @view p[:J]    
    @log H = J*ω #+ Hi    
    dx.θ = ω
    dx.ω = inv(J)*(x.T - ω×H) 
    return nothing
end

function fsw!(integrator)
    integrator.u.T = controller(integrator)
end

function controller(i)
    θ = @view i.u[:θ]
    ω = @view i.u[:ω]
    kp= @view i.p[:kp]
    kd= @view i.p[:kd]
    θr= @view i.p[:θr]
    ωr= @view i.p[:ωr]    

    T = kp*(θr-θ) + kd*(ωr-ω)    
    return T
end


x0 = ComponentArray(
    θ = zeros(3),
    ω = [pi/4,0,0], 
    T = zeros(3)   
    )

p = ComponentArray(    
    J = 1000.0*I(3),
    kp = 1000,
    kd = 1000,
    θr = zeros(3),
    ωr = zeros(3)
)
prob = ODEProblem(eom!,x0,(0,80),p,)
#=
saved_values = SavedValues(Float64,typeof(p))
save_func = (u,t,integrator) -> copy(integrator.p) # need to copy or all saved vals will be the end vals

save = SavingCallback(save_func,saved_values)
=#
fsw = PeriodicCallback(fsw!,.1)
sol = solve(prob,Tsit5(),callback = fsw)