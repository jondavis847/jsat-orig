function reactionWheelTorqueCommand!(S)

    # Negate body torque command since RWA applies torque internally
    S.u.fsw.output.rw.torque_command = S.p.fsw.output.rw.b_to_rw * -S.u.fsw.ac.ctrl.torque_command

    #Momentum Redistribtion and Zero Speed Avoidance
    momred!(S)

    RwCmdTqLim = S.u.fsw.output.rw.torque_command + S.u.fsw.output.rw.rwMomRdTq
    
    # Proportionally scale torque if torque command exceeds wheel capability    
    for idxWh = Base.OneTo(4)
        if abs(RwCmdTqLim[idxWh]) > S.p.fsw.output.rw.WhTqMax[idxWh]
            RwCmdTqLim = RwCmdTqLim / abs(RwCmdTqLim[idxWh]) * S.p.fsw.output.rw.WhTqMax[idxWh]
        end
    end
    S.u.fsw.output.rw.torque_command_lim = RwCmdTqLim
    S.u.fsw.output.rw.current_command  = RwCmdTqLim ./ S.p.actuators.rw.km
    #cmdCnt = uint16(CmdCurrent * 32768 / 7 + 32768))
    return nothing
end

function momred!(S)
    S.u.fsw.output.rw.prevMomAdj = S.u.fsw.output.rw.momAdj
    ## MOMRED is used for all wheel control modes to redistribute the wheel momentum in order to minimize the maximum momentum magnitude.
    # In addition, when the wheel speeds are close to zero, additional redistribution
    # is done to prevent the wheel speeds from lingering around zero.
    # (The wheel speeds do go through zero from time to time but they should
    # pass through reasonably quickly.) Keeping the speeds away from zero avoids
    # operating in the less accurate tachometer range
    # (where the lubricant distribution is uneven).
    
    # Make the null vector without changing the sign of the Wheel #1 momentum
    rwaAdjMom  = S.p.fsw.output.rw.nullVec .* S.u.actuators.rw.Hw
    
    # Find the momentum distribution for efficient momentum storage
    momMin  = minimum(rwaAdjMom)
    momMax  = maximum(rwaAdjMom)
    momAdj  = -0.5*(momMin + momMax)
    maxMomAdj = rwaAdjMom .+ momAdj
    
    #
    # ZERO SPEED AVOIDANCE
    # Since the computed wheel momenta after the first step may involve one or
    # more wheels running at speeds close to zero, further adjustment may be necessary.
    # Check for any speeds close to zero and calculate a new momentum distribution
    # using only POSITIVE (increasing momentum) adjustment for zero speed avoidance, if necessary.
    posMomAdj = wheelPosZeroAdj(maxMomAdj, S.u.actuators.rw.Hw, S.p.fsw.output.rw.rwaZeroMomThrd)
    # Starting from the resultant momentum distribution after the first step,
    # calculate a new momentum distribution using only NEGATIVE (decreasing momentum)
    # adjustment for zero speed avoidance, if necessary.
    negMomAdj = wheelNegZeroAdj(maxMomAdj, S.u.actuators.rw.Hw, S.p.fsw.output.rw.rwaZeroMomThrd)
    
    # Choose which (positive or negative) zero speed avoidance adjustment to apply.
    momAdj = (abs(posMomAdj) < abs(negMomAdj)) ? posMomAdj : negMomAdj
    
    S.u.fsw.output.rw.posMomAdj = posMomAdj
    S.u.fsw.output.rw.negMomAdj = negMomAdj
    # Avoid dithering in deadband
    S.u.fsw.output.rw.shiftDirCnt =  (momAdj * S.u.fsw.output.rw.prevMomAdj) < 0 ? S.u.fsw.output.rw.shiftDirCnt + 1 : 0
    
    
    if (( momAdj * S.u.fsw.output.rw.prevMomAdj) < 0) && (S.u.fsw.output.rw.shiftDirCnt <= S.p.fsw.output.rw.shiftDirThrd)
        momAdj = S.u.fsw.output.rw.prevMomAdj
    end
    
    #}
    
    # Add a ground commanded "null momentum" to avoid particular wheel speeds
    # that cause resonance in the telescope structure
    momAdj = momAdj + sign(momAdj)*S.p.fsw.output.rw.cmdNullMomZeroAvoid
    S.u.fsw.output.rw.momAdj = momAdj
    
    # Multiply by a gain to compute redistribution torque and limit it to a torque budget for momentum redistribution.
    rdTq = clamp.(S.p.fsw.output.rw.momRdGain * momAdj, -S.p.fsw.output.rw.whRdTqLim, S.p.fsw.output.rw.whRdTqLim)
    
    # Multiply by the actual null vector [1,-1,1,-1] to get redistribution wheel torques.
    # Note that, since these wheel torques have the ratio of the actual null vector,
    # the net torque on the spacecraft is theoretically zero;
    # however, a small net torque can actually exist due to wheel misalignment
    # and the variation of the actual wheel torque gains among the wheels
    S.u.fsw.output.rw.rwMomRdTq  = S.p.fsw.output.rw.nullVec * rdTq;
end# momred()

function wheelPosZeroAdj(PreMomAdj, RwMom, RwZeroMomThreshold)
    # Initialize zero speed avoidance momentum distribution to previously adjusted momentum distribution (4x1).
    Mom = PreMomAdj
    # Set AdjFound to TRUE; begin sweep for adjustment. AdjFound should become FALSE when all wheels are outside the zero speed zone defined by the threshold.
    AdjFound = true
    while AdjFound
        AdjFound = false
        for i = 1:4
            # Check whether wheel momenta are close to zero.
            if ( abs(Mom[i]) < RwZeroMomThreshold )
                # If any wheel momentum is found to be close to zero, AdjFound flag is set to TRUE so that there will be another sweep through all wheels to see if they all clear the zone of zero speed.
                AdjFound = true
                if (Mom[i] > 0)
                    # Mom[i] is positive and is close to zero. Increase all the momenta by the threshold value (null vector is [1,1,1,1]) so that Mom[i] becomes outside of the zone of zero speed.
                    for j = 1:4
                        Mom[j] = Mom[j] + RwZeroMomThreshold
                    end#for
                else
                    # Mom[i] is negative and is close to zero. Increase all the momenta by two times the threshold value (null vector is [1,1,1,1]) so that Mom[i] becomes outside of the zone of zero speed.
                    for j = 1:4
                        Mom[j] = Mom[j] + 2*RwZeroMomThreshold
                    end#for
                end#if
            end#if
        end#for
    end#while
    # Set output adjusted momentum (adjustment for all four wheels is the same, so the Wheel #1 value is returned).
    MomAdj = Mom[1] - RwMom[1]
end# wheelPosZeroAdj

function wheelNegZeroAdj(PreMomAdj, RwMom, RwZeroMomThreshold)
    # Initialize zero speed avoidance momentum distribution to previously adjusted momentum distribution (4x1).
    Mom = PreMomAdj
    # Set AdjFound to TRUE; begin sweep for adjustment. AdjFound should become FALSE when all wheels are outside the zero speed zone defined by the threshold.
    AdjFound = true
    while AdjFound
        AdjFound = false
        for i = 1:4
            # Check whether wheel momenta are close to zero.
            if abs(Mom[i]) < RwZeroMomThreshold
                # If any wheel momentum is found to be close to zero, AdjFound flag is set to TRUE so that there will be another sweep through all wheels to see if they all clear the zone of zero speed.
                AdjFound = true
                if Mom[i] > 0
                    # Mom[i] is positive and is close to zero. Decrease all the momenta by two times the threshold value (null vector is [1,1,1,1]) so that Mom[i] becomes outside of the zone of zero speed.
                    for j = 1:4
                        Mom[j] = Mom[j] - 2*RwZeroMomThreshold
                    end#for
                else
                    # Mom[i] is negative and is close to zero. Decrease all the momenta by the threshold value (null vector is [1,1,1,1]) so that Mom[i] becomes outside of the zone of zero speed.
                    for j = 1:4
                        Mom[j] = Mom[j] - RwZeroMomThreshold
                    end#for
                end#if
            end#if
        end#for
    end#while
    # Set output adjusted momentum (adjustment for all four wheels is the same, so the Wheel #1 value is returned).
    MomAdj = Mom[1] - RwMom[1]
end# wheelNegZeroAdj