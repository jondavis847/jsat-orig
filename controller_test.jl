using ModelingToolkit, IfElse
@variables t
function timer(func,dt;name)        
    @variable u(t) y(t)
    p = @parameters value = 0
    eq = [0 ~ y - value]
    ODESystem(eq, t; name=name, continuous_events = [(mod(t,dt) ~ 0) => (value = func())])
end




