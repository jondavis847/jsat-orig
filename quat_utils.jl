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

        function qmult