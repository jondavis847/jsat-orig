""" Actuators """

function actuators!(dx, x, p, t)
    reactionWheels!(dx, x, p, t)
    return nothing
end

function actuators_cb!(integrator)
    reactionWheels_cb!(integrator)
    mtb_cb!(integrator)
    return nothing
end

""" Reaction Wheels """

function reactionWheels!(dx, x, p, t)
    dx.actuators.rw.Hw = x.actuators.rw.Tw
    return nothing
end

function reactionWheels_cb!(S)
    S.u.actuators.rw.ω .= S.u.actuators.rw.Hw ./ S.p.actuators.rw.J
    S.u.actuators.rw.Hb .= hcat((S.u.actuators.rw.Hw .* eachcol(S.p.actuators.rw.a))...)

    #Back EMF
    S.u.actuators.rw.V_bemf .= S.p.actuators.rw.f_bemf .* abs.(S.u.actuators.rw.ω) 
    S.u.actuators.rw.V_motor .= S.u.fsw.output.rw.current_command .* S.p.actuators.rw.R_motor - S.u.actuators.rw.V_bemf
    S.u.actuators.rw.I_motor .= S.u.actuators.rw.V_motor ./ S.p.actuators.rw.R_motor
    S.u.actuators.rw.T_motor .= S.p.actuators.rw.km .* S.u.actuators.rw.I_motor

    #Torque curve applied to motor torque only?
    for i in eachindex(S.u.actuators.rw.Hw)
        if (S.u.actuators.rw.Hw[i] <= S.p.actuators.rw.knee[i])
            tqmax = -S.p.actuators.rw.H_m[i] * abs(S.u.actuators.rw.Hw[i]) + S.p.actuators.rw.H_b[i]
        else
            tqmax = -S.p.actuators.rw.T_m[i] * abs(S.u.actuators.rw.Hw[i]) + S.p.actuators.rw.T_b[i]
        end
        S.u.actuators.rw.T_motor[i] = clamp(S.u.actuators.rw.T_motor[i], -tqmax, tqmax)
    end

    #Ripple & Cogging torques    
    S.u.actuators.rw.T_cogging = sin.(S.u.actuators.rw.ω .* S.p.actuators.rw.n_poles ./ 2 .* S.t .+ S.p.actuators.rw.cogging_phase) * S.p.actuators.rw.cogging_amplitude
    S.u.actuators.rw.T_ripple = sin.(S.u.actuators.rw.ω .* S.p.actuators.rw.n_poles .* S.p.actuators.rw.n_phases .* S.t .+ S.p.actuators.rw.ripple_phase) * S.p.actuators.rw.ripple_coeff

    
    #friction
    S.u.actuators.rw.f_viscous = S.p.actuators.rw.f_viscous .* abs.(S.u.actuators.rw.ω)
    S.u.actuators.rw.f_windage = S.p.actuators.rw.f_windage .* S.u.actuators.rw.ω .^ 2
    for i in eachindex(S.u.actuators.rw.ω)
        if (abs(S.u.actuators.rw.ω[i]) < S.p.actuators.rw.stiction_speed[i])
            S.u.actuators.rw.f_stiction[i] = S.p.actuators.rw.f_stiction[i] .* tanh.(1 / S.p.actuators.rw.beta_s .* abs(S.u.actuators.rw.ω[i]))
            S.u.actuators.rw.T_friction[i] = sign(S.u.actuators.rw.ω[i]) * S.u.actuators.rw.f_stiction[i]
        else
            S.u.actuators.rw.f_stiction[i] = 0
            S.u.actuators.rw.T_friction[i] = sign(S.u.actuators.rw.ω[i]) * (S.p.actuators.rw.f_coulomb[i] + S.u.actuators.rw.f_viscous[i] + S.u.actuators.rw.f_windage[i])
        end
    end

    S.u.actuators.rw.Tw = S.u.actuators.rw.T_motor +
                          (S.u.actuators.rw.T_motor .* S.u.actuators.rw.T_ripple) +
                          S.u.actuators.rw.T_cogging -
                          S.u.actuators.rw.T_friction

    
    S.u.actuators.rw.Tb .= hcat((S.u.actuators.rw.Tw .* eachcol(S.p.actuators.rw.a))...)
end

""" Magnetic Torquer Bars """
function mtb_cb!(S)
    lut = linear_interpolation(S.p.actuators.mtb.current_lookup, S.p.actuators.mtb.moment_lookup) # stupid I have to do this every step instead of parameterizing
    S.u.actuators.mtb.Mm = lut(S.u.actuators.mtb.I) + S.p.actuators.mtb.residual_moment
    #S.u.actuators.mtb.Mm =  S.p.actuators.mtb.residual_moment
    S.u.actuators.mtb.Mb = S.p.actuators.mtb.mtb_to_brf * S.u.actuators.mtb.Mm
    S.u.actuators.mtb.Tb = cross(S.u.actuators.mtb.Mb, S.u.environments.geomagnetism.B_b)
    return nothing
end

