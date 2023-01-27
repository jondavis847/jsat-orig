includet("src\\jsat.jl")
includet("pace_parameters.jl")

"""Initial States"""

time = DateTime(2023,1,1,12,0,0)
epoch = Dates.datetime2julian(time)

x_orbit = ComponentArray(
    epoch =  epoch,
)

r0 =  [-6.32053381835774e6, -1.49222035669749e6, -2.77429961381375e6]
rmag0 = norm(r0)
v0 =  [2.56786681792351e3, 1.79527532306182e3, -6.82715713742553e3]

eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), epoch,eop_IAU1980)
r_ecef = eci_to_ecef*r0
lla = ecef_to_lla(r_ecef)

x_body = ComponentArray(    
    q = fakeNadir(r0,v0),
    ω = [0, -0.0011, 0], # [rad/sec],
    Hb = zeros(3),#body momentum
    Hi = zeros(3),#internal momentum
    Hs = zeros(3),#system momentum
    Te = zeros(3),
    Ti = zeros(3),
    eci_to_ecef = eci_to_ecef,
    r_eci = r0,
    r_ecef = r_ecef,
    lla = lla,
    v_eci = v0,
    rmag = rmag0
)
    
x_controller = ComponentArray(
    u = zeros(4),
    TqCmdBcs = zeros(3),    
    qr = fakeNadir(r0,v0),
    ωr = [0, -0.0011, 0], 
    attitudeError = zeros(3),
    rateError = zeros(3),
    integralError = zeros(3)
)

wheel_speeds = 100*ones(4)#zeros(4) #rad/sec
x_rw = ComponentArray(
    ω = SVector{4}(wheel_speeds), #wheel speed
    Tw = SVector{4}(zeros(4)), #wheel torque
    Tb = zeros(3,length(wheel_speeds)), #wheel torque in body frame
    Hw = p.actuators.rw.J .* wheel_speeds, # wheel momentum
    Hb = zeros(3,length(wheel_speeds)) # wheel momentum in body frame
)

x_mtb = ComponentArray(
    I = zeros(3),
    Mm = zeros(3),
    Mb = zeros(3),
    Tb = zeros(3)
)

x_actuators = ComponentArray(
    rw = x_rw,
    mtb = x_mtb
)

x_gravity = ComponentArray(a = zeros(3))
x_geomagnetism = ComponentArray(B_ecef = zeros(3),B_eci = zeros(3),B_b = zeros(3))
x_atmosphere = ComponentArray(ρ = 0.0)

x_environments = ComponentArray(
    geomagnetism = x_geomagnetism,
    gravity = x_gravity,
    atmosphere = x_atmosphere
)

x0 = ComponentArray(
    body = x_body,
    controller = x_controller,
    actuators = x_actuators,    
    orbit = x_orbit,
    environments = x_environments
)

#used for testing and benchmarking callback functions in REPL
S = ComponentArray(
    u = x0,
    p = p
)
""" ODE """
sol = simulate!(x0,p,(0,1000))