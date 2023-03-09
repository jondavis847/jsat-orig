includet("models\\dynamics.jl")
includet("models\\environments.jl")
includet("models\\orbit.jl")
includet("models\\actuators.jl")
includet("fsw\\fsw.jl")

using StaticArrays, DifferentialEquations, SatelliteToolbox

abstract type ModelValue end
Base.length(in::ModelValue) = length(in.value)
Base.one(in::ModelValue) = ones(size(in.value))

mutable struct NominalValue <: ModelValue
    value
end
mutable struct DispersedValue <: ModelValue
    value
    dist
end

Base.@kwdef struct InitialConditions
    Mode::ModelValue
    t::ModelValue
    r_eci::ModelValue
    v_eci::ModelValue
    q_i2b::ModelValue
    ω_b::ModelValue
    rwa_ω::ModelValue
end


""" Model """
function model!(dx, x, p, t)
    environments!(dx, x, p, t)
    actuators!(dx, x, p, t)
    eom!(dx, x, p, t)
end

function model_cb!(integrator)
    orbit_cb!(integrator)
    eom_cb!(integrator)
    environments_cb!(integrator)
    actuators_cb!(integrator)   # make sure act is after env (mtb) 
    fsw!(integrator)
end

lograte = 0.1
model_cb = PeriodicCallback(model_cb!, lograte, save_positions=(false, false))

function initModelValues(in, MC=false)
    T = Base.typename(typeof(in)).wrapper
    d = Dict() #initialize output
    v(x::NominalValue) = x.value
    v(x::DispersedValue) = !MC ? x.value : rand.(x.dist)    
    for k in keys(in)
        if in[k] isa T
            # if it's another CA, dig deeper
            merge!(d, Dict(k => initModelValues(in[k], MC)))
        elseif in[k] isa ModelValue
            merge!(d, Dict(k => v(in[k])))
            
        elseif in[k] isa Vector{Command} #reset commands when starting
            nt = [ComponentArray(f! = c.f!, t_start = c.t_start, duration = c.duration, occurring = false, occurred = false) for c in in[k]]
            merge!(d, Dict(k => nt))
        else
            #if it's something else, just put it back
            merge!(d, Dict(k => in[k]))
        end
    end
    return T(d)
end

function defineModel(ic)
    x_orbit = ComponentArray(
        epoch=ic.t,
        epoch_prev=Float64(0),
        sun_position_eci=SVector{3,Float64}(zeros(3)),
        sun_vector_eci=SVector{3,Float64}(zeros(3)),
        sun_vector_b=SVector{3,Float64}(zeros(3)),
        moon_position_eci=SVector{3,Float64}(zeros(3)),
        moon_vector_eci=SVector{3,Float64}(zeros(3)),
        moon_vector_b=SVector{3,Float64}(zeros(3)),
    )

    x_body = ComponentArray(
        q=ic.q_i2b,
        ω=ic.ω_b, # [rad/sec],
        Hb=SVector{3,Float64}(zeros(3)),#body momentum
        Hi=SVector{3,Float64}(zeros(3)),#internal momentum
        Hs=SVector{3,Float64}(zeros(3)),#system momentum
        Te=SVector{3,Float64}(zeros(3)),
        Ti=SVector{3,Float64}(zeros(3)),
        eci_to_ecef=SMatrix{3,3,Float64}(zeros(3, 3)),
        r_eci=ic.r_eci,
        r_ecef=SVector{3,Float64}(zeros(3)),
        lla=SVector{3,Float64}(zeros(3)),
        v_eci=ic.v_eci,
        v_ecef=SVector{3,Float64}(zeros(3)),
        rmag=Float64(0)
    )

    x_sun = ComponentArray(
        num_css_detect=Int64(0),
        sun_vector_desired=SVector{3,Float64}(zeros(3)),
        sun_vector_measured=SVector{3,Float64}(zeros(3)),
    )
    x_ad = ComponentArray(
        sun=x_sun,
        q_i2b=SVector{4,Float64}([0, 0, 0, 1]),
        ω_i2b=SVector{3,Float64}(zeros(3)),
    )

    x_ctrl = ComponentArray(
        attitude_error=SVector{3,Float64}(zeros(3)),
        attitude_error_lim=SVector{3,Float64}(zeros(3)),
        attitude_error_mag=Float64(0),
        integral_error=SVector{3,Float64}(zeros(3)),
        integral_error_lim=SVector{3,Float64}(zeros(3)),
        rate_error=SVector{3,Float64}(zeros(3)),
        rate_error_lim=SVector{3,Float64}(zeros(3)),
        torque_command=SVector{3,Float64}(zeros(3)),
        torque_command_b=SVector{3,Float64}(zeros(3)),
        torque_gyroscopic_ffwd=SVector{3,Float64}(zeros(3)),
        torque_oci_ffwd=SVector{3,Float64}(zeros(3)),
    )
    x_dyn = (
        Hs_b=SVector{3,Float64}(zeros(3)),
        Hs_b_mag=Float64(0),
        Hs_b_prev=SVector{3,Float64}(zeros(3)),
        ΔHs_b=SVector{3,Float64}(zeros(3)),
        ΔHs_b_mag=Float64(0),
    )

    x_logic = ComponentArray(
        gnc_mode=ic.Mode,
        gnc_mode_prev=ic.Mode,
    )

    x_tgt = ComponentArray(
        r_sc_to_site_eci=SVector{3,Float64}(zeros(3)),
        R_eci_to_brf=SMatrix{3,3,Float64}(zeros(3, 3)),
        q_ref=SVector{4,Float64}([0, 0, 0, 1]),
        q_nadir=SVector{4,Float64}([0, 0, 0, 1]),
        q_yaw_steering=SVector{4,Float64}([0, 0, 0, 1]),
        ω_ref=SVector{3,Float64}(zeros(3)),
        torque_oci_ffwd=SVector{3,Float64}(zeros(3))
    )

    x_ac = ComponentArray(
        ad=x_ad,
        ctrl=x_ctrl,
        dyn=x_dyn,
        logic=x_logic,
        tgt=x_tgt,
    )

    x_output_rw = ComponentArray(
        torque_command=SVector{4,Float64}(zeros(4)),
        torque_command_lim=SVector{4,Float64}(zeros(4)),
        current_command=SVector{4,Float64}(zeros(4)),
        rwMomRdTq=SVector{4,Float64}(zeros(4)),
        momAdj=Float64(0),
        prevMomAdj=Float64(0),
        posMomAdj=Float64(0),
        negMomAdj=Float64(0),
        shiftDirCnt=Int64(0),
    )
    x_output = ComponentArray(
        rw=x_output_rw,
    )

    x_commands = [ComponentArray(occurring=false,occurred=false)]

    x_fsw = ComponentArray(
        ac=x_ac,
        output=x_output,
        commands = x_commands,
    )

    x_rw = ComponentArray(
        ω=ic.rwa_ω, #wheel speed
        Tw=SVector{4,Float64}(zeros(4)), #wheel torque
        Tb=SMatrix{3,4,Float64}(zeros(3, 4)), #wheel torque in body frame
        Hw=SVector{4,Float64}(zeros(4)), # wheel momentum
        Hb=SMatrix{3,4,Float64}(zeros(3, 4)), # wheel momentum in body frame
        V_bemf=SVector{4,Float64}(zeros(4)),
        V_motor=SVector{4,Float64}(zeros(4)),
        I_motor=SVector{4,Float64}(zeros(4)),
        T_motor=SVector{4,Float64}(zeros(4)),
        T_cogging=SVector{4,Float64}(zeros(4)),
        T_ripple=SVector{4,Float64}(zeros(4)),
        f_viscous=SVector{4,Float64}(zeros(4)),
        f_windage=SVector{4,Float64}(zeros(4)),
        f_stiction=SVector{4,Float64}(zeros(4)),
        T_friction=SVector{4,Float64}(zeros(4)),
    )

    x_mtb = ComponentArray(
        I=SVector{3,Float64}(zeros(3)),
        Mm=SVector{3,Float64}(zeros(3)),
        Mb=SVector{3,Float64}(zeros(3)),
        Tb=SVector{3,Float64}(zeros(3))
    )

    x_actuators = ComponentArray(
        rw=x_rw,
        mtb=x_mtb
    )

    x_gravity = ComponentArray(
        a=zeros(3),
        gradient_torque_b=SVector{3,Float64}(zeros(3)),
    )
    x_geomagnetism = ComponentArray(
        B_ecef=SVector{3,Float64}(zeros(3)),
        B_eci=SVector{3,Float64}(zeros(3)),
        B_b=SVector{3,Float64}(zeros(3)))
    x_atmosphere = ComponentArray(ρ=Float64(0.0))

    x_environments = ComponentArray(
        geomagnetism=x_geomagnetism,
        gravity=x_gravity,
        atmosphere=x_atmosphere
    )

    ComponentArray(
        actuators=x_actuators,
        body=x_body,
        environments=x_environments,
        fsw=x_fsw,
        orbit=x_orbit,
    )
end
    #=
        rmag = norm(ic.r_eci)
        eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), epoch, p.environments.geomagnetism.eop_IAU1980)
        r_ecef = eci_to_ecef*ic.r_eci
        g = eci_to_ecef' * compute_g(p.environments.gravity.egm96_coefs, r_ecef, 16, 16)

        sunMOD = sun_position_i(epoch)
        mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),epoch)
        sun_position_eci = mod_to_j2000 * sunMOD
        sun_vector_eci = normalize(sun_position_eci - ic.r_eci)
        sun_vector_b = qvrot(ic.q_i2b,sun_vector_eci)

        moonMOD = moon_position_i(epoch)
        mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),epoch)
        moon_position_eci = mod_to_j2000 * moonMOD
        moon_vector_eci = normalize(moon_position_eci - ic.r_eci)
        moon_vector_b = qvrot(ic.q_i2b,moon_vector_eci)
    =#


function initModel(x, p)    
    x0 = initModelValues(x,p.config.montecarlo)
    p0 = initModelValues(p,p.config.montecarlo)

    #setup commands
    for i in eachindex(p0.fsw.commands)
        push!(x0.fsw.commands,ComponentArray(occurring = false,occurred = false))
    end

    #run once to setup all ic dependent values
    S = (u=x0, p=p0, t = 0, tprev = 0)
    model_cb!(S)
    return S
end