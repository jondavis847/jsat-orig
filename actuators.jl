""" Actuators """

function actuators!(dx, x, p, t)
    reactionWheels!(dx, x, p, t)
end

function actuators_cb!(integrator)
    reactionWheels_cb!(integrator)
    mtb_cb!(integrator)
end

""" Reaction Wheels """

function reactionWheels!(dx, x, p, t)    
    dx.actuators.rw.Hw = x.actuators.rw.Tw
end

function reactionWheels_cb!(S)    
    S.u.actuators.rw.Tw = S.p.actuators.rw.km .* S.u.controller.u
    S.u.actuators.rw.ω = S.u.actuators.rw.Hw ./ S.p.actuators.rw.J    
    S.u.actuators.rw.Tb = [(S.u.actuators.rw.Tw .* eachcol(S.p.actuators.rw.a))'...;]'
    S.u.actuators.rw.Hb = [(S.u.actuators.rw.Hw .* eachcol(S.p.actuators.rw.a))'...;]'
end

""" Magnetic Torquer Bars """
function mtb_cb!(S)
    lut = linear_interpolation(S.p.actuators.mtb.current_lookup,S.p.actuators.mtb.moment_lookup) # stupid I have to do this every step instead of parameterizing
    S.u.actuators.mtb.Mm = lut(S.u.actuators.mtb.I) + S.p.actuators.mtb.residual_moment
    #S.u.actuators.mtb.Mm =  S.p.actuators.mtb.residual_moment
    S.u.actuators.mtb.Mb = S.p.actuators.mtb.mtb_to_brf*S.u.actuators.mtb.Mm
    S.u.actuators.mtb.Tb = cross(S.u.actuators.mtb.Mb,S.u.environments.geomagnetism.B_b)
end

