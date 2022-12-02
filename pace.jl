includet("jsat.jl")

"""Define Parameters"""

mass = 1506.27
inertia = [1529.097 -58.0883 -26.71023
            -58.0883  1400.682 83.51491
            -26.71023  83.51491 2320.778]
cg =  [-1.311708,     -0.122079,    -0.0302493]

p_body = ComponentArray(
    J = inertia,
    invJ = inv(inertia),
    mass = mass,
    cg = cg,    
)

a1 =  [0.76604,0.64279,0]
a2 =  [0,0.64279,0.76604]
a3 =  [-0.76604,0.64279,0]
a4 =  [0,0.64279,-0.76604]
a = [a1,a2,a3,a4]

b_to_rw = [a'...;]

p_controller = ComponentArray(
    kp = [0.2853, 0.0288, 0.2853],
    kd = [0.7478, 0.2376, 0.7478],  
    ki = [0.0071, 0.0010, 0.0071],
    refrate = [0, -0.0011, 0],
    b_to_rw = b_to_rw,
    rw_to_b = b_to_rw',
    aeLim_MS    = [0.25, 0.75, 0.5] * pi/180,
    reLim_MS    = 1e6*ones(3),
    iaeLim_MS   = 0.0333*ones(3),
    EnableMinimax = 0,       # Flag to select distribution logic (0=DIS, 1=ENA)
    momRdGain   = 0.0,              # Momentum redistribution gain
    nullVec     = [1, -1, 1, -1],     # Null vector of mounting matrix, RwaToBcs
    whRdTqLim   = 0.225,            # Wheel torque limit for torque redistribution [Nm]
    cmdNullMomZeroAvoid = 0,        # Null momentum to avoid particular wheel speeds that cause resonance [Nm]
    rwaZeroMomThrd = 0,            # Wheel speed is considered close to zero if wheel momentum magnitude is below this value [Nms]
    shiftDirThrd = 1e3,             # Threshold for persistence counter value to exceed before actually changing momentum adjustment direction
    uijUD       = [    # upper diagonal terms of uij in Linf computation
        0.5077    0.0000   -0.5077   -0.5077   -1.0154   -0.5077
        -0.8520    0.0000    0.8520   -0.8520    0.0000   -0.8520
        0.5077    1.0154    0.5077    0.5077    0.0000   -0.5077
        ],        
    wijUD       = [# upper diagonal terms of wij in Linf computation
            0.3264    0.0000   -0.3264   -0.3264   -0.6527   -0.3264
           -0.3889    0.0000    0.3889   -0.3889    0.0000   -0.3889
            0.3264    0.6527    0.3264    0.3264    0.0000   -0.3264
            ],
        CoulombFricCoeff = zeros(4), # Coulomb friction coefficient of each wheel [Nm]
        ViscousFricCoeff = zeros(4), # Viscous friction coefficient of each wheel [Nms/rad]
        rwaDragGain      = 0,          # Scalar gain to compensate wheel drag
        EnableMomMax     = 0   # Flag to select enforcement of wheel momentum magnitude limit (0=DIS, 1=ENA)
  )

p_rw1 = ComponentArray(    
    km = 0.077, #motor constant
    J = 0.231, #wheel inertia
    a = a[1], #axis of rotation in reference frame 
    WhTqMax = 0.45,   # Max wheel torque (N*m)
    WhMomMax = 70,     # Max wheel momentum (N*m*s)   
)

p_rw2 = ComponentArray(    
    km = 0.077, #motor constant
    J = 0.231, #wheel inertia
    a = a[2], #axis of rotation in reference frame 
    WhTqMax = 0.45,   # Max wheel torque (N*m)
    WhMomMax = 70,     # Max wheel momentum (N*m*s)    
)

p_rw3 = ComponentArray(    
    km = 0.077, #motor constant
    J = 0.231, #wheel inertia
    a = a[3], #axis of rotation in reference frame    
    WhTqMax = 0.45,   # Max wheel torque (N*m)
    WhMomMax = 70,     # Max wheel momentum (N*m*s)  
)

p_rw4 = ComponentArray(    
    km = 0.077, #motor constant
    J = 0.231, #wheel inertia
    a = a[4], #axis of rotation in reference frame  
    WhTqMax = 0.45,   # Max wheel torque (N*m)
    WhMomMax = 70,     # Max wheel momentum (N*m*s)    
)

p_rw = [p_rw1,p_rw2,p_rw3,p_rw4]
     
p_gravity = ComponentArray(
    μ = 3.986004418e14,
    R = 6.378137e6,
    J2 =  1.08262668355e-3,
    J3 = -2.53265648533e-6,
    J4 = -1.61962159137e-6,
    J5 = -2.27296082869e-7,
    J6 =  5.40681239107e-7,
)

#p_geomagnetism = ComponentArray(
#    eop_IAU1980 = get_iers_eop()
#)

p_environments = ComponentArray(
    gravity = p_gravity,
    #geomagnetism = p_geomagnetism
)
p = ComponentArray(    
    body = p_body,
    controller = p_controller,    
    rw = p_rw,
    environments = p_environments,
)


"""Initial States"""

epoch = ComponentArray(
    Y = 2023.,
    M = 1.,
    D = 1.,
    h = 12.,
    m = 0.,
    s = 0.

)
x_orbit = ComponentArray(
    epoch =  epoch
)

r0 =  [-6.32053381835774e6, -1.49222035669749e6, -2.77429961381375e6]
rmag0 = norm(r0)
v0 =  [2.56786681792351e3, 1.79527532306182e3, -6.82715713742553e3]

x_body = ComponentArray(
    #θ = fakeNadir(r0,v0),
    q = fakeNadir(r0,v0),
    ω = [0, -0.0011, 0], # [rad/sec],
    H = zeros(3),
    Hi = zeros(3),
    Te = zeros(3),
    Ti = zeros(3),
    r = r0,
    v = v0,
    rmag = rmag0
)
    
x_controller = ComponentArray(
    u = zeros(4),
    TqCmdBcs = zeros(3),
    #θr = fakeNadir(r0,v0),
    qr = fakeNadir(r0,v0),
    ωr = [0, -0.0011, 0], 
    attitudeError = zeros(3),
    rateError = zeros(3),
    integralError = zeros(3)
)

ω = 0
x_rw1 = ComponentArray(
    ω = ω, #wheel speed
    Tw = 0, #wheel torque
    Tb = zeros(3), #wheel torque in body frame
    Hw = p.rw[1].J * ω, # wheel momentum
    Hb = zeros(3) # wheel momentum in body frame
)

ω = 0
x_rw2 = ComponentArray(
    ω = ω, #wheel speed
    Tw = 0, #wheel torque
    Tb = zeros(3), #wheel torque in body frame
    Hw = p.rw[2].J * ω, # wheel momentum
    Hb = zeros(3) # wheel momentum in body frame
)

ω = 0
x_rw3 = ComponentArray(
    ω = ω, #wheel speed
    Tw = 0, #wheel torque
    Tb = zeros(3), #wheel torque in body frame
    Hw = p.rw[3].J * ω, # wheel momentum
    Hb = zeros(3) # wheel momentum in body frame
)

ω = 0
x_rw4 = ComponentArray(
    ω = ω, #wheel speed
    Tw = 0, #wheel torque
    Tb = zeros(3), #wheel torque in body frame
    Hw = p.rw[4].J * ω, # wheel momentum
    Hb = zeros(3) # wheel momentum in body frame
)
x_rw = [x_rw1,x_rw2,x_rw3,x_rw4]

x_gravity = ComponentArray(
    a = zeros(3)
)

x_geomagnetism = ComponentArray(
    B = zeros(3)
)
x_environments = ComponentArray(
    geomagnetism = x_geomagnetism,
    gravity = x_gravity,
)

x0 = ComponentArray(
    body = x_body,
    controller = x_controller,
    rw = x_rw,    
    orbit = x_orbit,
    environments = x_environments
)


""" ODE """
sol = simulate!(x0,p,(0,5000))