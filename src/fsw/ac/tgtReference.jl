function tgtReference!(S)
    q_prev = S.u.fsw.ac.tgt.q_ref
    ω_prev = S.u.fsw.ac.tgt.ω_ref

    if (S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.MissionScience) ||
        (S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.DeltaV) ||
        (S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.DeltaH)

        S.u.fsw.ac.tgt.q_ref = S.u.fsw.ac.tgt.q_nadir        
        S.u.fsw.ac.tgt.ω_ref = S.p.fsw.ac.tgt.orbital_rate

    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.InertialReference

        S.u.fsw.ac.tgt.q_ref = S.p.fsw.ac.tgt.quat_train
        S.u.fsw.ac.tgt.ω_ref = @SVector zeros(3)

    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.LunarCal

        S.u.fsw.ac.tgt.q_ref = S.p.fsw.ac.tgt.quat_train
         
        q = S.u.fsw.ac.tgt.q_ref
        qdot = (q - q_prev)/0.1
        
        Χ = [q[4]*I(3)+scross(q[1:3]);-q[1:3]'] #fixme allocations
        ω = 2*Χ'*qdot
        S.u.fsw.ac.tgt.ω_ref = (ω + ω_prev)/2
    elseif S.u.fsw.ac.logic.gnc_mode == S.p.fsw.ac.logic.gnc_mode.SunSafe
        S.u.fsw.ac.ad.sun.sun_vector_desired = S.p.fsw.ac.ad.sun.sun_vector_desired #make ad.sun consistent, make this a quat? doesnt control roll axis
        S.u.fsw.ac.ad.sun.sun_vector_measured = S.u.orbit.sun_vector_b
        S.u.fsw.ac.tgt.q_ref = SVector{4,Float64}(0, 0, 0, 1)#comes from vector, not quat, maybe make this quat?
        S.u.fsw.ac.tgt.ω_ref = S.p.fsw.ac.ad.sun.rates_desired
    else
        S.u.fsw.ac.tgt.q_ref = SVector{4,Float64}(0, 0, 0, 1)
        S.u.fsw.ac.tgt.ω_ref = @SVector zeros(3)
    end
    if S.u.fsw.ac.tgt.q_ref'*q_prev < 0
        S.u.fsw.ac.tgt.q_ref = -S.u.fsw.ac.tgt.q_ref;
    end 

    S.u.fsw.ac.tgt.torque_oci_ffwd = @SVector zeros(3) # fix me at oci trq
    return nothing
end

scross(x) = SMatrix{3,3,Float64}(0, x[3], -x[2], -x[3], 0, x[1], x[2], -x[1], 0)

