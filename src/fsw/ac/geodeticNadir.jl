function geodeticNadir!(S)
    # find radius of curvature at site location
    # Rn = Rearth/sqrt(1-ecc^2*sin(lat)^2);
    # find position of site wrt to earth
    # r_EarthToSite_Ecef = Rn*[cos(lat)*cos(lon); cos(lat)*sin(lon); (1-ecc^2)*sin(lat)];
    # r_ScToSite_Eci = R_EcefToEci*(r_EarthToSite_Ecef - r_EarthToSc_Ecef);

    r_ScToSite_Eci = -S.u.body.eci_to_ecef' * S.u.body.lla[3] *
                     [cos(S.u.body.lla[1]) * cos(S.u.body.lla[2]),
                         cos(S.u.body.lla[1]) * sin(S.u.body.lla[2]),
                         sin(S.u.body.lla[1])]

    # construct basis vectors
    A3 = normalize(r_ScToSite_Eci)
    A2 = normalize(cross(-S.u.body.r_eci, S.u.body.v_eci))
    A1 = normalize(cross(A2, A3))
    # attitude matrix
    R_EciToBrf = [A1 A2 A3]'
    
    S.u.fsw.ac.tgt.r_sc_to_site_eci = r_ScToSite_Eci
    S.u.fsw.ac.tgt.R_eci_to_brf = R_EciToBrf
    S.u.fsw.ac.tgt.q_nadir= atoq(R_EciToBrf)
    return nothing
end
