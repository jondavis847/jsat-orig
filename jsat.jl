using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack  

everyStep(u,t,integrator) = mod(integrator.iter,1)
lograte = 0.1

function getsol(sol,names...)
    return map(sol.u) do u
        reduce(Base.maybeview,names; init=u)
    end
end

""" Model """
function model!(dx,x,p,t)
    actuators!(dx,x,p,t)
    eom!(dx,x,p,t)    
end

function gravity!(dx,x,p,t)
    @unpack μ,R,J2,J3,J4,J5,J6 = p.gravity
    @unpack x,y,z,r = x.body.position
    aj[1] = -(3/2)*J2*(μ/r^2)*(R/r)^2*[
        (1-5*(z/r)^2)*x/r
        (1-5*(z/r)^2)*y/r
        (3-5*(z/r)^2)*z/r
        ]

    aj[2] = -(1/2)*J3*(μ/r^2)*(R/r)^3*[
        5*(7*(z/r)^3 - 3*(z/r))*x/r
        5*(7*(z/r)^3 - 3*(z/r))*y/r
        3*(10*(z/r)^2 - 35/3*(z/r)^4 - 1)
        ]

    aj[3] = -(5/8)*J4*(μ/r^2)*(R/r)^4*[
        (3-42*(z/r)^2 + 63*(z/r)^4)*x/r
        (3-42*(z/r)^2 + 63*(z/r)^4)*y/r
        -(15-70*(z/r)^2 + 63*(z/r)^4)*z/r
        ]

    aj[4] = -(1/8)*J5*(μ/r^2)*(R/r)^5*[
        3*(35*(z/r)-210*(z/r)^3 + 231*(z/r)^5)*x/r
        3*(35*(z/r)-210*(z/r)^3 + 231*(z/r)^5)*y/r
        (15-315*(z/r)^2 + 945*(z/r)^4 - 693*(z/r)^6)
        ]

    aj[5] = (1/16)*J6*(μ/r^2)*(R/r)^6*[
        (35-945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6)*x/r
        (35-945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6)*y/r
        (2205*(z/r)^2 - 4851*(z/r)^4 + 3003*(z/r)^6 - 315)*z/r
        ]

    x.gravity.a = sum(aj)
end

gravity = PeriodicCallback(gravity!,lograte)

""" Body Translation """

function bodyTranslationCallback!(integrator)
    x.body.rmag = norm(x.body.r)
end

function bodyTranslation!(dx,x,p,t)
    @unpack r = dx.body
    @unpack rmag = dx.body

    dx.body.r = x.body.v
    dx.body.v = -p.gravity.μ/rmag*r + x.gravity.a 
end

""" Body Rotation """
function eom!(dx,x,p,t)
    bodyRotation!(dx,x,p,t)
end

function bodyRotationCallback!(integrator)        
    integrator.u.body.H = integrator.p.body.J*integrator.u.body.ω + integrator.u.body.Hi     
end
bodyRotationCallback = PeriodicCallback(bodyRotationCallback!,lograte)

function bodyRotation!(dx,x,p,t)  
    dx.body.θ = x.body.ω
    dx.body.ω = p.body.invJ*(x.body.Te - x.body.Ti - cross(x.body.ω,x.body.H))         
end

#function quaternionPropagation!(dx,x,p,t)
#    dx.body.q̇[1] = 
#end
""" Actuators """

### Reaction Wheels ###
function actuators!(dx,x,p,t)
    reactionWheels!(dx,x,p,t)
end

function reactionWheelCallback!(integrator)
    for i = 1:length(integrator.u.rw)
        integrator.u.rw[i].Tw = integrator.p.rw[i].km*integrator.u.controller.u[1]        
        integrator.u.rw[i].Hw = integrator.p.rw[i].J*integrator.u.rw[i].ω        
        integrator.u.rw[i].Tb = integrator.u.rw[i].Tw * integrator.p.rw[i].a
        integrator.u.rw[i].Hb = integrator.u.rw[i].Hw * integrator.p.rw[i].a
    end
    integrator.u.body.Ti = sum([integrator.u.rw[i].Tb for i in 1:length(integrator.u.rw)])
    integrator.u.body.Hi = sum([integrator.u.rw[i].Hb for i in 1:length(integrator.u.rw)])
end

rw = PeriodicCallback(reactionWheelCallback!,lograte)

function reactionWheels!(dx,x,p,t)
    for i = 1:length(x.rw)
        dx.rw[i].ω = x.rw[i].Tw/p.rw[i].J        
    end
end

""" FSW """

function fsw!(integrator)
    attitudeError!(integrator)
    controller!(integrator)
end
fswrate = 0.1
fsw = PeriodicCallback(fsw!,fswrate)

function attitudeError!(integrator)
    #integrator.u.controller.θr = norm(-integrator.u.body.r)
    integrator.u.controller.θr = [0,0,0]
    integrator.u.controller.ωr = [0,0,0]
end

function controller!(integrator)
    integrator.u.controller.u = (
        integrator.p.controller.kp.*(integrator.u.controller.θr-integrator.u.body.θ)
         + integrator.p.controller.kd.*(integrator.u.controller.ωr-integrator.u.body.ω))
    for i = 1:3
        if integrator.u.controller.u[i] > 7
            integrator.u.controller.u[i] = 7
        end
        if integrator.u.controller.u[i] < -7
            integrator.u.controller.u[i] = -7
        end
    end
end

