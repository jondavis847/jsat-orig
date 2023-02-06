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
    acMassMomCalc!(integrator)
    selectLogic!(integrator)
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