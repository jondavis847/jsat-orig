function acMassMomCalc!(S)
    S.u.fsw.ac.dyn.Hs_b_prev = S.u.fsw.ac.dyn.Hs_b
    S.u.fsw.ac.dyn.Hs_b = S.u.body.Hb + S.u.body.Hi #fix me sensor inputs, not truth
    S.u.fsw.ac.dyn.Hs_b_mag = norm(S.u.fsw.ac.dyn.Hs_b)
    
    S.u.fsw.ac.dyn.ΔHs_b = S.u.fsw.ac.dyn.Hs_b - S.u.fsw.ac.dyn.Hs_b_prev
    S.u.fsw.ac.dyn.ΔHs_b_mag = norm(S.u.fsw.ac.dyn.ΔHs_b)
    return nothing
end