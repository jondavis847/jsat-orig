abstract type Actuators end
abstract type Bodies end
abstract type Controllers end
abstract type Dynamics end
abstract type Frame end
abstract type Parameters end
abstract type Sensors end
abstract type States end

struct SpacecraftStates <: States
    actuators::Vector{Actuators}
    bodies::Vector{Bodies}
    controllers::Vector{Controllers}
    sensors::Vector{Sensors}
    dynamics::Vector{Dynamics}
end

struct Spacecraft
    actuators::Vector{Actuators}
    bodies::Vector{Bodies}
    controllers::Vector{Controllers}
    sensors::Vector{Sensors}
    dynamics::Vector{Dynamics}
end

function spacecraft_to_component_array(sc)
end

function 

function eom!(dx,x,t,p)
    #dx.sensors = sensors(dx,x,t,p)
    #dx.controllers = controllers(dx,x,t,p)
    #dx.actuators = actuators(dx,x,t,p)
    dx.dynamics = dynamics(dx,x,t,p)
end

