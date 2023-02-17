function attitudeController!(S)
    # Attitude controller to stabilize 3-axis of the spacecraft                       
    # Control error

    S.u.fsw.ac.ctrl.attitude_error = qtov(qmult(qinv(S.u.fsw.ac.tgt.q_ref), S.u.fsw.ac.ad.q_i2b))#S.u.fsw.ac.ad.q_i2b))
    if norm(S.u.fsw.ac.ctrl.attitude_error) > pi # rotate in minimum-course direction
        S.u.fsw.ac.ctrl.attitude_error = qtov(qmult(qinv(-S.u.fsw.ac.tgt.q_ref), S.u.fsw.ac.ad.q_i2b))#S.u.fsw.ac.ad.q_i2b))
    end
    #S.u.fsw.ac.ctrl.rate_error = S.u.fsw.ac.ad.ω_i2b - S.u.fsw.ac.tgt.ω_ref
    S.u.fsw.ac.ctrl.rate_error = S.u.fsw.ac.ad.ω_i2b - S.u.fsw.ac.tgt.ω_ref
    S.u.fsw.ac.ctrl.integral_error = S.u.fsw.ac.ctrl.integral_error + S.u.fsw.ac.ctrl.attitude_error * 0.1 # hardcoded for 10Hz processing

    # set control attitude error and integral of control attitude
    # error to zero if target or estimated quaternions are not valid
    # attitude quaternions
    #tlmpkt.attCtrlValidFlag = ValidInvalid.Valid;
    #tol = this.unitQuat_tol;
    #if abs(vecnorm(refQuat)-1)>tol || abs(vecnorm(estQuat)-1)>tol
    #   ctrlAttErr = zeros(3,1);
    #   ctrlIntAttErr_ = zeros(3,1);
    #   tlmpkt.attCtrlValidFlag = ValidInvalid.Invalid;
    #end

    # GNC mode dependent controller logic
    if S.u.fsw.ac.logic.gnc_mode == p.fsw.ac.logic.gnc_mode.RateNull

        S.u.fsw.ac.ctrl.attitude_error_lim = S.u.fsw.ac.ctrl.attitude_error  # not using it. Don't care.
        S.u.fsw.ac.ctrl.rate_error_lim = S.u.fsw.ac.ctrl.rate_error # NO limiter
        S.u.fsw.ac.ctrl.integral_error_lim = zeros(3)  # not using it. Don't care.

        # No gyroscopic, no OCI ffwd torque
        S.u.fsw.ac.ctrl.torque_command_b = S.p.fsw.ac.ctrl.moi * (-p.fsw.ac.ctrl.Kd_RN .* S.u.fsw.ac.ctrl.rate_error_lim)
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = zeros(3)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = zeros(3)

    elseif S.u.fsw.ac.logic.gnc_mode == p.fsw.ac.logic.gnc_mode.SunSafe

        # proportional control error
        if (S.u.fsw.ac.ad.sun.num_css_detect == S.p.fsw.ac.ad.sun.css_detect.many_css) ||
           (S.u.fsw.ac.ad.sun.num_css_detect == S.p.fsw.ac.ad.sun.css_detect.two_css) ||
           (S.u.fsw.ac.ad.sun.num_css_detect == S.p.fsw.ac.ad.sun.css_detect.one_css)
            # if the measured sun vector is valid, compute the attitude
            # error from it and the desired sun vector
            attErrDir = normalize(cross(S.u.orbit.sun_vector_b, S.u.fsw.ac.ad.sun.sun_vector_desired))
            S.u.fsw.ac.ctrl.attitude_error = acos(dot(
                normalize(S.u.orbit.sun_vector_b),
                normalize(S.u.fsw.ac.ad.sun.sun_vector_desired)
            )) * attErrDir
        else
            # Otherwise, if its invalid (most often due to being in
            # eclipse, but possibly due to other reasons), override disable
            # the attitude portion of the controller via zeroing the error
            S.u.fsw.ac.ctrl.attitude_error = zeros(3)
        end
        S.u.fsw.ac.ctrl.attitude_error_lim = S.u.fsw.ac.ctrl.attitude_error

        # derivative control error
        S.u.fsw.ac.ctrl.rate_error = S.u.fsw.ac.ad.ω_i2b - S.u.fsw.ac.tgt.ω_ref 
        S.u.fsw.ac.ctrl.rate_error_lim = S.u.fsw.ac.ctrl.rate_error

        S.u.fsw.ac.ctrl.integral_error = zeros(3)
        S.u.fsw.ac.ctrl.integral_error_lim = zeros(3) # not using it. Don't care.

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = - (S.p.fsw.ac.ctrl.moi * (
            S.p.fsw.ac.ctrl.Kp_SS .* S.u.fsw.ac.ctrl.attitude_error +
            S.p.fsw.ac.ctrl.Kd_SS .* S.u.fsw.ac.ctrl.rate_error
        ))
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = zeros(3)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = zeros(3)

    elseif S.u.fsw.ac.logic.gnc_mode == p.fsw.ac.logic.gnc_mode.MissionScience

        S.u.fsw.ac.ctrl.attitude_error_lim = clamp.(S.u.fsw.ac.ctrl.attitude_error, -S.p.fsw.ac.ctrl.aeLim_MS, S.p.fsw.ac.ctrl.aeLim_MS)
        S.u.fsw.ac.ctrl.rate_error_lim = clamp.(S.u.fsw.ac.ctrl.rate_error, -S.p.fsw.ac.ctrl.reLim_MS, S.p.fsw.ac.ctrl.reLim_MS)

        # reset integrated error at mode entry
        if S.u.fsw.ac.logic.gnc_mode_prev != S.u.fsw.ac.logic.gnc_mode
            S.u.fsw.ac.ctrl.integral_error = @SVector zeros(3)
        end

        # switch between PD/PID ctrl
        if any((abs.(S.u.fsw.ac.ctrl.attitude_error) .> S.p.fsw.ac.ctrl.aeLim_MS) .| (abs.(S.u.fsw.ac.ctrl.rate_error) .> S.p.fsw.ac.ctrl.reLim_MS))
            # PD ctrl
            S.u.fsw.ac.ctrl.integral_error_lim = @SVector zeros(3)
        else
            # PID ctrl
            S.u.fsw.ac.ctrl.integral_error_lim = clamp.(S.u.fsw.ac.ctrl.integral_error, -S.p.fsw.ac.ctrl.iaeLim_MS, S.p.fsw.ac.ctrl.iaeLim_MS)
        end
    
        ang_accel = -(S.p.fsw.ac.ctrl.Kp_MS .* S.u.fsw.ac.ctrl.attitude_error_lim
                    +S.p.fsw.ac.ctrl.Kd_MS .* S.u.fsw.ac.ctrl.rate_error_lim
                    +S.p.fsw.ac.ctrl.Ki_MS .* S.u.fsw.ac.ctrl.integral_error_lim)

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = S.p.fsw.ac.ctrl.moi * ang_accel
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = -cross(S.u.fsw.ac.ad.ω_i2b,S.u.fsw.ac.dyn.Hs_b)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = S.u.fsw.ac.tgt.torque_oci_ffwd

    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.InertialReference
        # Compute error quaternion                  
        eigvec = normalize(S.u.fsw.ac.ctrl.attitude_error)
        ctrlEulerAng = norm(S.u.fsw.ac.ctrl.attitude_error)
        # euler angle control limit
        ctrlEulerAngLim = clamp(ctrlEulerAng, -S.p.fsw.ac.ctrl.aeAngLim_nonlin_IR, S.p.fsw.ac.ctrl.aeAngLim_nonlin_IR)

        # Assign Values
        S.u.fsw.ac.ctrl.attitude_error_lim = eigvec * sin(ctrlEulerAngLim / 2) # scaled error quaternion                 
        S.u.fsw.ac.ctrl.rate_error_lim = clamp.(S.u.fsw.ac.ctrl.rate_error, -S.p.fsw.ac.ctrl.reLim_IR, S.p.fsw.ac.ctrl.reLim_IR)

        # reset integrated error at mode entry
        if S.u.fsw.ac.logic.gnc_mode != S.u.fsw.ac.logic.gnc_mode_prev
            S.u.fsw.ac.ctrl.integral_error = zeros(3)
        end

        # switch between PD/PID ctrl
        if any(
            (abs.(S.u.fsw.ac.ctrl.attitude_error) .> S.p.fsw.ac.ctrl.aeLim_IR) ||
            (abs.(S.u.fsw.ac.ctrl.rate_error) .> S.p.fsw.ac.ctrl.reLim_IR)
        )
            # PD ctrl
            S.u.fsw.ac.ctrl.integral_error_lim = zeros(3)
        else
            # PID ctrl
            S.u.fsw.ac.actrl.integral_error_lim = clamp.(S.u.fsw.integral_error, -S.p.fsw.ac.ctrl.iaeLim_IR, S.p.fsw.ac.ctrl.iaeLim_IR)
        end
        ang_accel = -(S.p.fsw.ac.ctrl.Kp_IR .* S.u.fsw.ac.ctrl.attitude_error_lim
        + S.p.fsw.ac.ctrl.Kd_IR .* S.u.fsw.ac.ctrl.rate_error_lim
        + S.p.fsw.ac.ctrl.Ki_IR .* S.u.fsw.ac.ctrl.integral_error_lim)

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = mean(diag(S.p.fsw.ac.ctrl.moi)) .* ang_accel
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = zeros(3)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = S.u.fsw.ac.tgt.torque_oci_ffwd

    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.LunarCal

        S.u.fsw.ac.ctrl.attitude_error_lim = clamp.(S.u.fsw.ac.ctrl.attitude_error, -S.p.fsw.ac.ctrl.aeLim_LC, S.p.fsw.ac.ctrl.aeLim_LC)
        S.u.fsw.ac.ctrl.rate_error_lim = clamp.(S.u.fsw.ac.ctrl.rate_error, -S.p.fsw.ac.ctrl.reLim_LC, S.p.fsw.ac.ctrl.reLim_LC)
        S.u.fsw.ac.ctrl.integral_error_lim = zeros(3)# Don't care. Not used.

        # PDD design but no acceleration term.                    
        ang_accel = -(S.p.fsw.ac.ctrl.Ki_LC .* S.p.fsw.ac.ctrl.attitude_error_lim
        +S.p.fsw.ac.ctrl.Kp_LC .* S.p.fsw.ac.ctrl.rate_error_lim)

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = S.p.fsw.ac.ctrl.moi .* angAccel
        # Add gyroscopic feed-forward torque
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = -cross(S.u.fsw.ac.ad.ω_i2b, S.u.fsw.ac.dyn.Hs_b)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = zeros(3)

    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.DeltaV

        S.u.fsw.ac.ctrl.attitude_error_lim = S.u.fsw.ac.ctrl.attitude_error
        S.u.fsw.ac.ctrl.rate_error_lim = S.u.fsw.ac.ctrl.rate_error

        # reset integrated error at mode entry
        if S.u.fsw.ac.logic.gnc_mode != S.u.fsw.ac.logic.gnc_mode_prev
            S.u.fsw.ac.ctrl.integral_error = zeros(3)
        end

        S.u.fsw.ac.ctrl.integral_error_lim = S.u.fsw.ac.ctrl.integral_error

        # PID control
        ang_accel = -(S.p.fsw.ac.ctrl.Kp_DV .* S.u.fsw.ac.ctrl.attitude_error_lim
        + S.p.fsw.ac.ctrl.Kd_DV .* S.u.fsw.ac.ctrl.rate_error_lim
        + S.p.fsw.ac.ctrl.Ki_DV .* S.u.fsw.ac.ctrl.integral_error_lim)

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = S.p.fsw.ac.ctrl.moi .* angAccel
        # Add gyroscopic feed-forward torque
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = -cross(S.u.fsw.ac.ad.ω_i2b, S.u.fsw.ac.dyn.Hs_b)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = zeros(3)


    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.DeltaH

        S.u.fsw.ac.ctrl.attitude_error_lim = S.u.fsw.ac.ctrl.attitude_error
        S.u.fsw.ac.ctrl.rate_error_lim = S.u.fsw.ac.ctrl.rate_error
        S.u.fsw.ac.ctrl.integral_error_lim = zeros(3)

        # PD control
        ang_accel = -(S.p.fsw.ac.ctrl.Kp_DH .* S.u.fsw.ac.ctrl.attitude_error_lim
        + S.p.fsw.ac.ctrl.Kd_DH .* S.u.fsw.ac.ctrl.rate_error_lim)

        # TqCmdBcs computation
        S.u.fsw.ac.ctrl.torque_command_b = S.p.fsw.ac.ctrl.moi .* angAccel
        # Add gyroscopic feed-forward torque
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd = -cross(S.u.fsw.ac.ad.ω_i2b, S.u.fsw.ac.dyn.Hs_b)
        S.u.fsw.ac.ctrl.torque_oci_ffwd = zeros(3)

    end
    # Add gyroscopic & OCI feed-fodard torque
    S.u.fsw.ac.ctrl.torque_command =
        S.u.fsw.ac.ctrl.torque_command_b +
        S.u.fsw.ac.ctrl.torque_gyroscopic_ffwd +
        S.u.fsw.ac.ctrl.torque_oci_ffwd

    S.u.fsw.ac.ctrl.attitude_error_mag = norm(S.u.fsw.ac.ctrl.attitude_error)
    S.u.fsw.ac.dyn.Hs_b_mag = norm(S.u.fsw.ac.dyn.Hs_b)
    S.u.fsw.ac.dyn.ΔHs_b = S.u.fsw.ac.dyn.Hs_b - S.u.fsw.ac.dyn.Hs_b_prev
    S.u.fsw.ac.dyn.ΔHs_b_mag = norm(S.u.fsw.ac.dyn.ΔHs_b)

    # pop last measurement and add new value
    #this.EstSysMomBcs_MagPrev = [this.EstSysMomBcs_MagPrev(2:end); EstSysMomBcs_mag_temp];
    # compute meant
    #tlmpkt.EstSysMomBcs_AvgMag = single(mean(this.EstSysMomBcs_MagPrev));
    # Update discrete states
    return nothing
end

