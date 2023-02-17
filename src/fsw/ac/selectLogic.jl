function selectLogic!(S)
    #placeholder
    q = ( S.u.body.q[4] >=0 ) ? S.u.body.q : -S.u.body.q
    S.u.fsw.ac.ad.q_i2b = q
    S.u.fsw.ac.ad.ω_i2b = S.u.body.ω
    return nothing
end