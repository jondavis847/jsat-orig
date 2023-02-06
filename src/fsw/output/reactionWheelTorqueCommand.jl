function reactionWheelTorqueCommand!(S)

    # Negate body torque command since RWA applies torque internally
    S.u.fsw.output.rw.torque_command = S.p.fsw.ac.ctrl.b_to_rw * -S.u.fsw.ac.ctrl.torque_command

    # Proportionally scale torque if torque command exceeds wheel capability
    RwCmdTqLim = S.u.fsw.output.rw.torque_command
    for idxWh = Base.OneTo(4)
        if abs(RwCmdTqLim[idxWh]) > S.p.actuators.rw.WhTqMax[idxWh]
            RwCmdTqLim = RwCmdTqLim / abs(RwCmdTqLim[idxWh]) * S.p.actuators.rw.WhTqMax[idxWh]
        end
    end
    S.u.fsw.output.rw.torque_command_lim = RwCmdTqLim
    S.u.fsw.output.rw.current_command  = RwCmdTqLim ./ p.actuators.rw.km
    #cmdCnt = uint16(CmdCurrent * 32768 / 7 + 32768))
    return nothing
end