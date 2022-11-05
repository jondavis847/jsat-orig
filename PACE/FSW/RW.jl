module RW_M
include("Parameters.jl")

function wheelCommand(wheelSpeeds, brfTorqueCommand)

    # Negate body torque command since RWA applies torque internally
    wheelTorqueCommands = -brfTorqueCommand
    # Compute wheel momenta in wheel hyperframe
    wheelMomentums = FSW.RW.J .* wheelSpeeds

    # use minimum norm pseudo inverse to map body frame torque
    # command into wheel frame
    wheelTorqueCommand = FSW.RW.BcsToRwa * wheelTorqueCommands

    # Zero-speed Avoidance & momentum REDistribution (ZARED)
    wheelMomDistTorque = zeros(4) #this.momred(wheelMomentums);

    # Feed-forward compensate wheel drag
    wheelDragCompTorque = zeros(4) #this.rwaDrag(wheelSpeeds);

    # Synthesize torque due to Controller & ZARED & Drag
    wheelTotalTorqueCommand = wheelTorqueCommand + wheelMomDistTorque + wheelDragCompTorque

    # Proportionally scale torque if torque command exceeds wheel capability
    wheelLimitedTorqueCommand = wheelTotalTorqueCommand
    for wheel = 1:4
        if abs(wheelLimitedTorqueCommand[wheel]) > FSW.RW.MaxTorque[wheel]
            wheelLimitedTorqueCommand = wheelLimitedTorqueCommand / abs(wheelLimitedTorqueCommand[wheel]) * FSW.RW.MaxTorque[wheel]
        end
    end

    # Zero the torque command if any wheel is saturated and torque
    # command will not reduce saturation
    Saturated = false
    for wheel = 1:4
        #   check if Mom exceeds thresh && check if torque will not improve the situation
        if (abs(wheelMomentums[wheel]) >= FSW.RW.MaxMomentum[wheel]) && (wheelMomentums[wheel] * FSW.RW.MaxTorque[wheel] > 0)
            Saturated = true
        end
    end
    if Saturated
        wheelLimitedTorqueCommand = zeros(4)
    end
    # Converting the desired torque into a motor current
    # Converting the motor current to an integer number of counts
    wheelCurrentCommand = wheelLimitedTorqueCommand ./ FSW.RW.KMotor
    #cmdCnt = uint16(wheelCurrentCommand * 32768 / 7 + 32768);
    return wheelCurrentCommand #cmdCnt
end

end