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

function yawSteering!(S)
    ## Yaw steering angle calculation
    w   = -7.29211585494e-5 # Earth rotation (spin) rate, [rad/s]
    P   = S.u.body.r_ecef
    V   = S.u.body.v_ecef
    lat,lon,alt = S.u.body.lla
    # Nadir vector in ECEF
    N   = -[cos(lat)*cos(lon)
           cos(lat)*sin(lon)
           sin(lat)]
    # Effect of the Earth rotation on velScEcef --> velScEci
    V_  = V + w*cross(P,SA[0, 0, 1])
    
    # Yaw rotation angle
    # tan ? = N?(V×V’) / ((V?V’) – (N?V)(N?V’))
    num = dot(N, cross(V, V_))
    den = dot(V,V_) - dot(N,V)*dot(N,V_)
    ysAng = atan(num, den)

    S.u.fsw.ac.tgt.q_yaw_steering= SA[(0.5*sin(ysAng)*SA[0,0,1])..., 0.5*cos(ysAng)]
    return nothing
end
#=
function refQuatRateProc(S)
    # REFQUATRATEPROC generates mode-dependent reference quaternion 
    # (inertial to body) and angular velocity, given ndrQuat, biasQuat,
    # quatTrain, and GNCmode.
    
    persistent refQPrev wPrev
    if isempty(refQPrev); refQPrev = [0 0 0 0]'; end
    if isempty(wPrev); wPrev = [0 0 0]'; end
    
    ## Reference quaternion/rate
    switch gncMode
        case {GNC_Mode.MissionScience,GNC_Mode.DeltaH,GNC_Mode.DeltaV}
            refQ = qmult(gdNdrQuat,biasQuat);
            if refQ'*refQPrev<0
                refQ = -refQ;
            end 
            w0 = [0 -AC_Target_TBL.orbitalRateMag 0]'; # [rad/sec]
            refW  = qvrot(w0,biasQuat);
        case GNC_Mode.InertialRef
            refQ = quatTrain;
            if refQ'*refQPrev<0
                refQ = -refQ;
            end 
            refW = zeros(3,1);
        case GNC_Mode.LunarCal
            # w0 = [0 0.0021816616 0]'; # 450 arcsec/s = 0.00218166156499 [rad/sec]
            # qdot = qmult(qinv(refQPrev),refQuat)/0.1;
            refQ = quatTrain;
            if refQ'*refQPrev<0
                refQ = -refQ;
            end 
            qdot = (refQ - refQPrev)/0.1;
            q = refQ;
            Xi = [q(4)*eye(3)+scross(q(1:3));-q(1:3)'];
            w = 2*Xi'*qdot;
            refW = (w + wPrev)/2;
        otherwise
            refQ = [0 0 0 1]';
            refW = zeros(3,1);
    end
    refQPrev = refQ;
    wPrev = refW;
=#