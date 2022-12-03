""" Actuators """

function actuators!(dx, x, p, t)
    reactionWheels!(dx, x, p, t)
end

function actuators_cb!(integrator)
    reactionWheels_cb!(integrator)
end

### Reaction Wheels ###

function reactionWheels!(dx, x, p, t)    
    dx.rw.ω = x.rw.Tw ./ p.rw.J    
end

function reactionWheels_cb!(S)
    
    S.u.rw.Tw = S.p.rw.km .* S.u.controller.u
    S.u.rw.Hw = S.p.rw.J .* S.u.rw.ω
    S.u.rw.Tb = [(S.u.rw.Tw .* eachcol(S.p.rw.a))'...;]
    S.u.rw.Hb = [(S.u.rw.Hw .* eachcol(S.p.rw.a))'...;]
    
    S.u.body.Ti = sum(S.u.rw.Tb,dims=2)
    S.u.body.Hi = sum(S.u.rw.Hb,dims=2)
end

