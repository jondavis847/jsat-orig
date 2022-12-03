""" FSW """
fswrate = 0.1

function fsw!(integrator)
    attitudeError!(integrator)
    controller!(integrator)
    reactionWheelTorqueCommand!(integrator)
end

fsw = PeriodicCallback(fsw!, fswrate)

function attitudeError!(S)
    #S.u.controller.θr = fakeNadir(S.u.body.r,S.u.body.v)
    S.u.controller.qr = fakeNadir(S.u.body.r, S.u.body.v)
    S.u.controller.ωr = S.p.controller.refrate

    S.u.controller.attitudeError = qtov(qmult(qinv(S.u.controller.qr), S.u.body.q))
    if norm(S.u.controller.attitudeError) > pi # rotate in minimum-course direction
        S.u.controller.attitudeError = qtov(qmult(qinv(-S.u.controller.qr), S.u.body.q))
    end
    S.u.controller.rateError = S.u.body.ω - S.u.controller.ωr
    S.u.controller.integralError = S.u.controller.integralError + S.u.controller.attitudeError * 0.1 # hardcoded for 10Hz processing
    #S.u.controller.attitudeError = S.u.controller.θr-S.u.body.θ
    #S.u.controller.rateError = S.u.controller.ωr-S.u.body.ω
end

function controller!(S)
    # Assign values               
    S.u.controller.attitudeError = limitLowerUpper(S.u.controller.attitudeError, -S.p.controller.aeLim_MS, S.p.controller.aeLim_MS)
    S.u.controller.rateError = limitLowerUpper(S.u.controller.rateError, -S.p.controller.reLim_MS, S.p.controller.reLim_MS)

    # reset integrated error at mode entry
    #if LogicTlm.Mode ~= this.lastMode
    #    ctrlIntAttErr_     = zeros(3,1);
    #end

    # switch between PD/PID ctrl
    if any((abs.(S.u.controller.attitudeError) .> S.p.controller.aeLim_MS) .| (abs.(S.u.controller.rateError) .> S.p.controller.reLim_MS))
        # PD ctrl
        S.u.controller.integralError = @SVector zeros(3)
    else
        # PID ctrl
        S.u.controller.integralError = limitLowerUpper(S.u.controller.integralError, -S.p.controller.iaeLim_MS, S.p.controller.iaeLim_MS)
    end

    # TqCmdBcs computation
    tqCmdBcs = S.p.body.J * (
        -S.p.controller.kp .* S.u.controller.attitudeError 
        -S.p.controller.kd .* S.u.controller.rateError
        -S.p.controller.ki .* S.u.controller.integralError)


    # Add gyroscopic & OCI feed-forward torque
    S.u.controller.TqCmdBcs = tqCmdBcs - cross(S.u.body.ω, S.u.body.H)# + tgt.OciFfwdTorque;

end

function reactionWheelTorqueCommand!(S)

    # Negate body torque command since RWA applies torque internally
    rwaTqCmd_Brf = -S.u.controller.TqCmdBcs

    rwTqCmd_w = S.p.controller.b_to_rw * rwaTqCmd_Brf

    # Proportionally scale torque if torque command exceeds wheel capability
    RwCmdTqLim = rwTqCmd_w
    for idxWh = 1:4
        if abs(RwCmdTqLim[idxWh]) > S.p.rw.WhTqMax[idxWh]
            RwCmdTqLim = RwCmdTqLim / abs(RwCmdTqLim[idxWh]) * S.p.rw.WhTqMax[idxWh]
        end
    end

    S.u.controller.u  = RwCmdTqLim ./ p.rw.km
    #cmdCnt = uint16(CmdCurrent * 32768 / 7 + 32768))
end

function fakeNadir(r, v)
    #A3 = -normalize(S.u.body.r)
    #A2 = normalize(cross(S.u.body.r,-S.u.body.v))
    #A1 = cross(a2,a3)
    #R_EciToBrf = [A1 A2 A3]';

    #Attempting not to allocate?
    R_EciToBrf = [cross(normalize(cross(r, -v)), -normalize(r)) normalize(cross(r, -v)) -normalize(r)]
    #q_EciToGdrf = atoq(R_EciToBrf)
    #angles = qtoe(q_EciToGdrf)

    R = RotMatrix{3,Float64}(R_EciToBrf)
    #a = RotZYX(R)
    #angles = Rotations.params(a)
    quat = QuatRotation(R)
    q = Rotations.params(quat)
    return @SVector [q[2],q[3],q[4],q[1]] #make it scalar last
end