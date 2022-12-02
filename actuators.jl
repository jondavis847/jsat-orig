""" Actuators """

function actuators!(dx, x, p, t)
    reactionWheels!(dx, x, p, t)
end

function actuators_cb!(integrator)
    reactionWheels_cb!(integrator)
end

### Reaction Wheels ###

function reactionWheels!(dx, x, p, t)
    for i = 1:length(x.rw)
        dx.rw[i].ω = x.rw[i].Tw / p.rw[i].J
    end
end

function reactionWheels_cb!(S)
    for i = 1:length(S.u.rw)
        S.u.rw[i].Tw = S.p.rw[i].km * S.u.controller.u[i]
        S.u.rw[i].Hw = S.p.rw[i].J * S.u.rw[i].ω
        S.u.rw[i].Tb = S.u.rw[i].Tw * S.p.rw[i].a
        S.u.rw[i].Hb = S.u.rw[i].Hw * S.p.rw[i].a
    end
    S.u.body.Ti = sum([S.u.rw[i].Tb for i in 1:length(S.u.rw)])
    S.u.body.Hi = sum([S.u.rw[i].Hb for i in 1:length(S.u.rw)])
end

