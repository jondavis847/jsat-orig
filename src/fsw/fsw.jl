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

#mutable struct Command
struct Command
    f!
    t_start::Float64 #seconds
    duration::Float64 #seconds
    #occurred::Bool
    #occurring::Bool
    #Command(f,t_start,duration) = new(f,t_start,duration,false,false)
end

Base.length(::Command) = 1

function commands!(S)    
    for i in eachindex(S.p.fsw.commands)
        if (S.t >= S.p.fsw.commands[i].t_start) && !S.u.fsw.commands[i].occurred
            if S.u.fsw.commands[i].occurring && (S.t >= (S.p.fsw.commands[i].t_start + S.p.fsw.commands[i].duration))                
                S.u.fsw.commands[i].occurred = true                
                S.u.fsw.commands[i].occurring = false
                continue
            end
            S.u.fsw.commands[i].occurring = true
            S.p.fsw.commands.f!(S)
        end                
    end
end