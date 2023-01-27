using LinearAlgebra
#inverse quaternion,  btime 132.775 ns (3 allocations: 288 bytes)
qinv(q) = q .* SA[-1, -1, -1, 1] / norm(q)^2

#quaternion to rotation vector,  #btime 197.586 ns (4 allocations: 320 bytes)
qtov(q) = 2 * atan(norm(q[1:3]), q[4]) * normalize(q[1:3])

#rotation matrix to quaternion
function atoq(A) # btime 548.634 ns (19 allocations: 928 bytes)   

    e1 = @SVector [0, 1.0 + A[1, 1] - A[2, 2] - A[3, 3]]
    e2 = @SVector [0, 1.0 - A[1, 1] + A[2, 2] - A[3, 3]]
    e3 = @SVector [0, 1.0 - A[1, 1] - A[2, 2] + A[3, 3]]
    e4 = @SVector [0, 1.0 + A[1, 1] + A[2, 2] + A[3, 3]]
    tmp = @SVector [sqrt(maximum(e1)),
        sqrt(maximum(e2)),
        sqrt(maximum(e3)),
        sqrt(maximum(e4))]

    f1 = abs.(tmp)
    qmax = maximum(f1)
    if f1 == qmax
        q1 = 0.5 * tmp[1]
        rq = 0.25 / q1
        q = @SVector [
            q1,
            rq * (A[1, 2] + A[2, 1]),
            rq * (A[1, 3] + A[3, 1]),
            rq * (A[2, 3] - A[3, 2]),
        ]
    elseif f2 == qmax
        q2 = 0.5 * tmp[2]
        rq = 0.25 / q2
        q = @SVector [
            rq * (A[1, 2] + A[2, 1]),
            q2,
            rq * (A[2, 3] + A[3, 2]),
            rq * (A[3, 1] - A[1, 3])
        ]
    elseif f3 == qmax
        q3 = 0.5 * tmp[3]
        rq = 0.25 / q3
        q = @SVector [
            rq * (A[3, 1] + A[1, 3]),
            rq * (A[3, 2] + A[2, 3]),
            q3,
            rq * (A[1, 2] - A[2, 1])
        ]
    else
        q4 = 0.5 * tmp[4]
        rq = 0.25 / q4
        q = @SVector [
            rq * (A[2, 3] - A[3, 2]),
            rq * (A[3, 1] - A[1, 3]),
            rq * (A[1, 2] - A[2, 1]),
            q4
        ]
    end

    # Normalize quaternion output
    return normalize(q)
end

function qtoe(q)
    # returns the euler angles (in radians) corresponding to a quaternions
    w = q[4]
    x = q[1]
    y = q[2]
    z = q[3]

    ysqr = y * y

    t0 = 2.0 * (w * x + y * z)
    t1 = 1.0 - 2.0 * (x * x + ysqr)
    x_ang = atan(t0, t1)

    t2 = 2.0 * (w * y - z * x)
    if (t2 > 1.0)
        t2 = 1.0
    end
    if (t2 < -1.0)
        t2 = -1.0
    end
    y_ang = asin(t2)

    t3 = 2.0 * (w * z + x * y)
    t4 = 1.0 - 2.0 * (ysqr + z * z)
    z_ang = atan(t3, t4)

    return @SVector [x_ang, y_ang, z_ang]

end

function qmult(q_A2B, q_B2C) # 162.533 ns (3 allocations: 288 bytes)
    # q_A2C = qmult(q_A2B,q_B2C, varargin)
    #
    # Quaternion multiplication:   q_A2B * q_B2C    = q_A2C 
    #                           A(q_B2C) * A(q_A2B) = A(q_A2C)
    #
    # this multiplication implementation follows the \odot convention from 
    # Markley, F. Landis, and John L. Crassidis. "Fundamentals of Spacecraft 
    # Attitude Determination and Control." (2014): 978-1.  page 37
    #
    # let quaternions be denoted q = [qv; qs] and p = [pv; ps]
    # then   q \odot p = [ ps*qv + qs*pv + cross(qv,pv);...
    #                              qs*ps - dot(qv,pv)  ];
    #

    q_A2C = @SVector [
        q_B2C[4] * q_A2B[1] + q_B2C[3] * q_A2B[2] - q_B2C[2] * q_A2B[3] + q_B2C[1] * q_A2B[4],
        -q_B2C[3] * q_A2B[1] + q_B2C[4] * q_A2B[2] + q_B2C[1] * q_A2B[3] + q_B2C[2] * q_A2B[4],
        q_B2C[2] * q_A2B[1] - q_B2C[1] * q_A2B[2] + q_B2C[4] * q_A2B[3] + q_B2C[3] * q_A2B[4],
        -q_B2C[1] * q_A2B[1] - q_B2C[2] * q_A2B[2] - q_B2C[3] * q_A2B[3] + q_B2C[4] * q_A2B[4]
    ]


    q_A2C = normalize(q_A2C)

    if q_A2C[4] < 0
        q_A2C = -q_A2C
    end

    return q_A2C
end # qmult()

function qvrot(q,v,Transform = true)
    #can be transform or rotation
    #transform: physical orientation is constant, but represented in different frames
    #rotation: physical orientation is changing, represented in the same frame
    vmag = norm(v)
    u_v = normalize(v)    
    v_aug = append!(u_v,0)
    if Transform
        q = qinv(q)
    end

    q = qmult(q,qmult(q,v_aug))    
    return vmag*view(q,1:3)
end

function getsol(sol, names...)
    return map(sol.u) do u
        reduce(Base.maybeview, names; init=u)
    end
end

function limitLowerUpper!(val, lower, upper)
    for i in eachindex(val)
        if val[i] > upper[i]
            val[i] = upper[i]
        elseif val[i] < lower[i]
            val[i] = lower[i]
        end

        #val[i] = val[i] > upper[i] ? upper[i] : val[i]
        #val[i] = val[i] < lower[i] ? lower[i] : val[i]
    end
    return val
end

