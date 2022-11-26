using ComponentArrays, DifferentialEquations, LinearAlgebra, UnPack, Rotations

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

function eom!(dx,x,p,t)
    bodyTranslation!(dx,x,p,t)
    bodyRotation!(dx,x,p,t)
end


function gravity!(integrator)
    @unpack μ,R,J2,J3,J4,J5,J6 = integrator.p.gravity
    x,y,z = @views integrator.u.body[:r]
    r = @views integrator.u.body[:rmag]

    aj2 = -(3/2)*J2*(μ/r^2)*(R/r)^2*[
        (1-5*(z/r)^2)*x/r
        (1-5*(z/r)^2)*y/r
        (3-5*(z/r)^2)*z/r
        ]

    aj3 = -(1/2)*J3*(μ/r^2)*(R/r)^3*[
        5*(7*(z/r)^3 - 3*(z/r))*x/r
        5*(7*(z/r)^3 - 3*(z/r))*y/r
        3*(10*(z/r)^2 - 35/3*(z/r)^4 - 1)
        ]

    aj4 = -(5/8)*J4*(μ/r^2)*(R/r)^4*[
        (3-42*(z/r)^2 + 63*(z/r)^4)*x/r
        (3-42*(z/r)^2 + 63*(z/r)^4)*y/r
        -(15-70*(z/r)^2 + 63*(z/r)^4)*z/r
        ]

    aj5 = -(1/8)*J5*(μ/r^2)*(R/r)^5*[
        3*(35*(z/r)-210*(z/r)^3 + 231*(z/r)^5)*x/r
        3*(35*(z/r)-210*(z/r)^3 + 231*(z/r)^5)*y/r
        (15-315*(z/r)^2 + 945*(z/r)^4 - 693*(z/r)^6)
        ]

    aj6 = (1/16)*J6*(μ/r^2)*(R/r)^6*[
        (35-945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6)*x/r
        (35-945*(z/r)^2 + 3465*(z/r)^4 - 3003*(z/r)^6)*y/r
        (2205*(z/r)^2 - 4851*(z/r)^4 + 3003*(z/r)^6 - 315)*z/r
        ]

    integrator.u.gravity.a = aj2+aj3+aj4+aj5+aj6
end

gravity = PeriodicCallback(gravity!,lograte)


function atoq(A)
#=
% This function computes the attitude quaternion corresponding
% to a given attitude matrix...
%
%  INPUT
%         A   Input attitude matrix.
%  OUTPUT
%         q   Output attitude quaternion.
%
%  EXTERNAL REFERENCES:  NONE.
%
%function [q]=atoq(A);

% Check norm of input Direction Cosine Matrix A
% if abs(abs(det(A))-1) > 1e-12
%     error('LIB:atoq','Input attitude matrix is not an orthonormal matrix.');
% end

% Now we can proceed ... =#
    q = zeros(4);
    e1=[0;1.0+A[1,1]-A[2,2]-A[3,3]]
    e2=[0;1.0-A[1,1]+A[2,2]-A[3,3]]
    e3=[0;1.0-A[1,1]-A[2,2]+A[3,3]]
    e4=[0;1+A[1,1]+A[2,2]+A[3,3]]
    q[1]=sqrt(maximum(e1))
    q[2]=sqrt(maximum(e2))
    q[3]=sqrt(maximum(e3))
    q[4]=sqrt(maximum(e4))
    f1=[abs(q[1]),abs(q[2]),abs(q[3]),abs(q[4])]
    qmax=maximum(f1)
    if abs(q[1]) == qmax
        q[1]=0.5*q[1]
        rq=0.25/q[1]
        q[2]=rq*(A[1,2]+A[2,1])
        q[3]=rq*(A[1,3]+A[3,1])
        q[4]=rq*(A[2,3]-A[3,2])
    elseif abs(q[2]) == qmax
        q[2]=0.5*q[2]
        rq=0.25/q[2]
        q[1]=rq*(A[1,2]+A[2,1])
        q[3]=rq*(A[2,3]+A[3,2])
        q[4]=rq*(A[1,3]-A[1,3])
    elseif abs(q[3]) == qmax
        q[3]=0.5*q[3]
        rq=0.25/q[3]
        q[1]=rq*(A[1,3]+A[1,3])
        q[2]=rq*(A[3,2]+A[2,3])
        q[4]=rq*(A[1,2]-A[2,1])
    else
        q[4]=0.5*q[4]
        rq=0.25/q[4]
        q[1]=rq*(A[2,3]-A[3,2])
        q[2]=rq*(A[1,3]-A[1,3])
        q[3]=rq*(A[1,2]-A[2,1])
    end
        
    # Normalize quaternion outpu
    return normalize(q)
end
    
function qtoe(q)
    # returns the euler angles (in radians) corresponding to a quaternions
        w = q[4];
        x = q[1];
        y = q[2];
        z = q[3];
        
        ysqr = y * y;
        
        t0 = 2.0 * (w * x + y * z);
        t1 = 1.0 - 2.0 * (x * x + ysqr);
        x_ang = atan(t0, t1);
        
        t2 = 2.0 * (w * y - z * x);
        if(t2 > 1.0)
            t2 = 1.0;  
        end
        if(t2 < -1.0)
            t2 = -1.0;
        end
        y_ang = asin(t2);
        
        t3 = 2.0 * (w * z + x * y);
        t4 = 1.0 - 2.0 * (ysqr + z * z);
        z_ang = atan(t3, t4);

        return [x_ang, y_ang, z_ang]
        
    end

function fakeNadir(r,v)
    #A3 = -normalize(integrator.u.body.r)
    #A2 = normalize(cross(integrator.u.body.r,-integrator.u.body.v))
    #A1 = cross(a2,a3)
    #R_EciToBrf = [A1 A2 A3]';
    
    #Attempting not to allocate?
    R_EciToBrf = [cross(normalize(cross(r,-v)),-normalize(r))  normalize(cross(r,-v))  -normalize(r) ]
    #q_EciToGdrf = atoq(R_EciToBrf)
    #angles = qtoe(q_EciToGdrf)

    R = RotMatrix{3,Float64}(R_EciToBrf)
    a = RotZYX(R)
    angles = Rotations.params(a)    

    return angles
end

""" Body Translation """

function bodyTranslation!(dx,x,p,t)
    dx.body.r = x.body.v
    dx.body.v = -p.gravity.μ/(norm(x.body.r)^3)*x.body.r + x.gravity.a 
end

""" Body Rotation """

function bodyRotationCallback!(integrator)        
    #integrator.u.body.H = integrator.p.body.J*integrator.u.body.ω + integrator.u.body.Hi     
    integrator.u.body.ω = integrator.p.body.invJ*(integrator.u.body.H -integrator.u.body.Hi)
end
bodyRotationCallback = PeriodicCallback(bodyRotationCallback!,lograte)

function bodyRotation!(dx,x,p,t)  
    dx.body.θ = x.body.ω
    dx.body.H = x.body.Te - x.body.Ti - cross(x.body.ω,x.body.H)
    #dx.body.ω = p.body.invJ*(x.body.Te - x.body.Ti - cross(x.body.ω,x.body.H))         
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
    individual_commands = integrator.p.controller.b_to_rw * integrator.u.controller.u
    for i = 1:length(integrator.u.rw)
        integrator.u.rw[i].Tw = integrator.p.rw[i].km*individual_commands[i]        
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
    integrator.u.controller.θr = fakeNadir(integrator.u.body.r,integrator.u.body.v)
    integrator.u.controller.ωr = integrator.p.controller.refrate
    integrator.u.controller.attitudeError = integrator.u.controller.θr-integrator.u.body.θ
    integrator.u.controller.rateError = integrator.u.controller.ωr-integrator.u.body.ω
end

function controller!(integrator)
    integrator.u.controller.u = -(
        integrator.p.controller.kp.*(integrator.u.controller.attitudeError)
         + integrator.p.controller.kd.*(integrator.u.controller.rateError))
    for i = 1:3
        if integrator.u.controller.u[i] > 7
            integrator.u.controller.u[i] = 7
        end
        if integrator.u.controller.u[i] < -7
            integrator.u.controller.u[i] = -7
        end
    end
end

