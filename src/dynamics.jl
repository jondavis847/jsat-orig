""" True Dynamics """
function eom!(dx, x, p, t)
    bodyTranslation!(dx, x, p, t)
    bodyRotation!(dx, x, p, t)
    return nothing
end

function eom_cb!(integrator)
    bodyRotation_cb!(integrator)
    bodyTranslation_cb!(integrator)
    return nothing
end

function time_cb!(S)
    S.u.orbit.epoch_prev = S.u.orbit.epoch
    S.u.orbit.epoch = S.u.orbit.epoch + (S.t - S.tprev) / 86400.0
    return nothing
end

# Body Translation

function bodyTranslation!(dx, x, p, t)
    # No longer propagating this way, using ObritPropagator from satellite SatelliteToolbox
    # This way produces incorrect orbit, raising altitude by 300km in 3 days.
    #dx.body.r_eci = x.body.v_eci    
    #dx.body.v_eci = x.environments.gravity.a
end

function bodyTranslation_cb!(S)

    orbp = init_orbit_propagator(Val(:J2), OrbitStateVector(S.u.orbit.epoch_prev,
        SVector{3,Float64}(S.u.body.r_eci),
        SVector{3,Float64}(S.u.body.v_eci),
        SVector{3,Float64}(S.u.environments.gravity.a)))

    r, v = propagate!(orbp, S.t - S.tprev)
    S.u.body.r_eci = r
    S.u.body.v_eci = v

    S.u.body.eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), S.u.orbit.epoch, S.p.environments.geomagnetism.eop_IAU1980)
    S.u.body.r_ecef = S.u.body.eci_to_ecef * S.u.body.r_eci
    S.u.body.v_ecef = S.u.body.eci_to_ecef * S.u.body.v_eci
    S.u.body.rmag = norm(S.u.body.r_eci)
    #S.u.body.lla = ecef_to_lla(S.u.body.r_ecef)
    S.u.body.lla .= ecef_to_geodetic(S.u.body.r_ecef)
    return nothing
end

# Body Rotation 

function bodyRotation!(dx, x, p, t)
    q = x.body.q
    # Crassidis/Markley, Eq 3.20 wrt Eq 2.88
    Q = @SMatrix [
        q[4] -q[3] q[2]
        q[3] q[4] -q[1]
        -q[2] q[1] q[4]
        -q[1] -q[2] -q[3]
    ]

    
    dx.body.q = 0.5 * Q * x.body.ω
    dx.body.Hb = x.body.Te - x.body.Ti - cross(x.body.ω, x.body.Hs)
    return nothing
end

function bodyRotation_cb!(S)
    S.u.body.Hs = S.u.body.Hb + S.u.body.Hi
    S.u.body.ω = inv(S.p.body.J) * S.u.body.Hb
    if S.u.body.q[4] < 0
        S.u.body.q = -S.u.body.q
    end

    S.u.body.Te = S.u.actuators.mtb.Tb + S.u.environments.gravity.gradient_torque_b
    S.u.body.Ti = sum(S.u.actuators.rw.Tb, dims=2)
    S.u.body.Hi = sum(S.u.actuators.rw.Hb, dims=2)
    return nothing
end
