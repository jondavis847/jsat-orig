includet("ac\\acMassMomCalc.jl")
includet("ac\\attitudeController.jl")
includet("ac\\geodeticNadir.jl")
includet("ac\\tgtReference.jl")
includet("ac\\selectLogic.jl")
includet("output\\mtbTorqueCommand.jl")
includet("output\\reactionWheelTorqueCommand.jl")

""" FSW """

function fsw!(integrator)            
    ac!(integrator)
    output!(integrator)
    return nothing
end

function ac!(integrator)
    commands!(integrator)
    acMassMomCalc!(integrator)
    selectLogic!(integrator)
    yawSteering!(integrator)
    geodeticNadir!(integrator)
    tgtReference!(integrator)
    attitudeController!(integrator)
    return nothing
end

function output!(integrator) #formally known as do, but conflicts with julia do loops
    reactionWheelTorqueCommand!(integrator)
    mtbTorqueCommand!(integrator)
    return nothing
end

mutable struct Command
    f!
    t_start::Float64 #seconds
    duration::Float64 #seconds
    occurred::Bool
    occurring::Bool
    Command(f,t_start,duration) = new(f,t_start,duration,false,false)
end

Base.length(::Command) = 1

function commands!(S)
    for c in S.p.fsw.commands            
        if (S.t >= c.t_start) && !c.occurred                  
            if c.occurring && (S.t >= (c.t_start + c.duration))                
                c.occurred = true                
                c.occurring = false
                continue
            end
            c.occurring = true
            c.f!(S)
        end                
    end
end