using ModelingToolkit, DifferentialEquations, Plots, IfElse
@variables t 

function step(time;name)    
    x = @variables u(t)=0 [output = true,irreducible = true]    
    p = @parameters time=time height=0        
    eq = [0 ~ u - height]
    ODESystem(eq, t; name=name, continuous_events = [(t - time ~ 0) => (height ~ 1)])
end

function steps(times;name)        
    n = length(times)
    x = @variables u(t)[1:n]=zeros(n) [output = true,irreducible = true]
    p = @parameters times[1:n]=times heights[1:n]=zeros(n)        

    eq = [0 ~ u[i] - heights[i] for i in 1:n]
        
    cond = [t - times[i] ~ 0 for i in 1:n]
    affect = [heights[i] ~ 1 for i in 1:n]
    ODESystem(eq,t,[x...;],[p...;]; name=name, continuous_events = cond => affect)
end

function step_ifelse(time;name)
    x = @variables u(t)=0 [output = true,irreducible = true]    
    p = @parameters time=time
    eq = [0 ~ u - IfElse.ifelse(t > time, 1., 0. )]
    ODESystem(eq, t; name=name)
end

function steps_ifelse(times;name)
    n = length(times)
    @variables u(t)[1:n]=zeros(n) [output = true,irreducible = true]    
    @parameters times[1:n]=times
    eq = [u[i] ~ IfElse.ifelse( times[i] < t, 1., 0. ) for i in 1:n]
    ODESystem(eq, t; name=name)
end