using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations, StaticArrays, Dates, PlotlyJS
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

includet("model.jl")
includet("utils\\plotting.jl")
includet("utils\\quaternion.jl")
includet("utils\\simulation.jl")