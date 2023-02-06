using UnPack, SatelliteToolbox
includet("orbit.jl")

"""Environments """

function environments_cb!(S)
    #gravity_cb!(integrator)
    #geomagnetism_cb!(integrator)
    stb_gravity_cb!(S)
    if S.p.config.geomagnetism
        stb_geomagnetism_cb!(S)
    end
    if S.p.config.atmosphere
        stb_atmosphere_cb!(S)
    end
    return nothing
end

""" Gravity """
function stb_gravity_cb!(S)
    S.u.environments.gravity.a = SVector{3}(S.u.body.eci_to_ecef' * compute_g(S.p.environments.gravity.egm96_coefs, S.u.body.r_ecef, 16, 16))
    return nothing
end
#over load length of gravity model to support component arrays
Base.length(in::GravityModel_Coefs{Float64}) = 0

""" Geomagnetism """
function stb_geomagnetism_cb!(S)
    #IGRF is currently 3 times slower than dipole (3.5us vs 9.4ns, 26 vs 12 allocations)
    decimal_year = decyear(Dates.julian2datetime(S.u.orbit.epoch))

    B_ned = igrf(decimal_year, S.u.body.lla[3], S.u.body.lla[1], S.u.body.lla[2], Val(:geodetic))
    S.u.environments.geomagnetism.B_ecef = 1e-9 * ned_to_ecef(B_ned, S.u.body.lla...)

    #becef =  SVector{3}(geomag_dipole(x0.body.r_ecef,decimal_year))
    #S.u.environments.geomagnetism.B_ecef = SVector{3}(1e-9*becef)
    #S.u.environments.geomagnetism.B_ecef = igrf13syn(1,decimal_year,1,S.u.body.lla[3]/1000, 180-S.u.body.lla[1]*180/pi, S.u.body.lla[2]*180/pi) * 1e-9

    S.u.environments.geomagnetism.B_eci = SVector{3}(S.u.body.eci_to_ecef' * S.u.environments.geomagnetism.B_ecef)
    S.u.environments.geomagnetism.B_b = qvrot(qinv(S.u.body.q), S.u.environments.geomagnetism.B_eci)
    return nothing
end

#need to overload EOP length to work with component arrays
Base.length(in::EOPData_IAU1980) = 0

""" Atmospheric Drag """
function stb_atmosphere_cb!(S)
    S.u.environments.atmosphere.ρ = expatmosphere(S.u.body.lla[3])
    return nothing
end