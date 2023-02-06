function  mtbTorqueCommand!(S)
    MomUnloadMagFld = cross(S.u.body.Hs, S.u.environments.geomagnetism.B_b)
    Dipole_Bcs = S.p.actuators.mtb.dipole_gain .* MomUnloadMagFld
    Dipole_Mtb = S.p.actuators.mtb.mtb_to_brf' * Dipole_Bcs + S.p.actuators.mtb.dipole_bias

   # Re-scaling to compensate for the saturated MTB axis:
    max_dipole = maximum(abs.(Dipole_Mtb)) 
    if max_dipole > S.p.actuators.mtb.dipole_maxlinear
        ScaleFactor = S.p.actuators.mtb.dipole_maxlinear / max_dipole
    else
        ScaleFactor = 1.0
    end
    
    Dipole_Mtb  = Dipole_Mtb * ScaleFactor
    S.u.actuators.mtb.I = S.p.actuators.mtb.Moment2Current_Slope .* Dipole_Mtb + S.p.actuators.mtb.Moment2Current_Intercept    
    return nothing
end