includet("src\\jsat.jl")
includet("pace_parameters.jl")

"""Initial States"""

r0 =  [-6.32053381835774e6, -1.49222035669749e6, -2.77429961381375e6]
rmag0 = norm(r0)
v0 =  [2.56786681792351e3, 1.79527532306182e3, -6.82715713742553e3]

time = DateTime(2023,1,1,12,0,0)
epoch = Dates.datetime2julian(time)
eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), epoch,eop_IAU1980)
r_ecef = eci_to_ecef*r0
g = eci_to_ecef' * compute_g(p.environments.gravity.egm96_coefs, r_ecef, 16, 16)
q = [-0.021179442940975148, 0.5496477692711645, 0.16029818602936372, 0.8195994463685526]

sunMOD = sun_position_i(epoch)
mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),epoch)
sun_position_eci = mod_to_j2000 * sunMOD
sun_vector_eci = normalize(sun_position_eci - r0)
sun_vector_b = qvrot(q,sun_vector_eci)

moonMOD = moon_position_i(epoch)
mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),epoch)
moon_position_eci = mod_to_j2000 * moonMOD
moon_vector_eci = normalize(moon_position_eci - r0)
moon_vector_b = qvrot(q,moon_vector_eci)


x_orbit = ComponentArray(
    epoch =  Float64(epoch),
    epoch_prev = Float64(epoch-0.1/86400.0),
    sun_position_eci = SVector{3,Float64}(sun_position_eci),
    sun_vector_eci = SVector{3,Float64}(sun_vector_eci),
    sun_vector_b = SVector{3,Float64}(sun_vector_b),
    moon_position_eci = SVector{3,Float64}(moon_position_eci),
    moon_vector_eci = SVector{3,Float64}(moon_vector_eci),
    moon_vector_b = SVector{3,Float64}(moon_vector_b),
)

#lla = Vector{Float64}(undef,3)
lla = collect(ecef_to_geodetic(r_ecef)) #returns a tuple, so need to preallocate and broadcast to vector

J = [1529.097 -58.0883 -26.71023; -58.0883 1400.682 83.51491; -26.71023 83.51491 2320.778]
ω = [0, -0.0011, 0]
Hb = J*ω

x_body = ComponentArray(    
    q = SVector{4,Float64}(q),
    ω = SVector{3,Float64}(ω), # [rad/sec],
    Hb = SVector{3,Float64}(Hb),#body momentum
    Hi = SVector{3,Float64}(zeros(3)),#internal momentum
    Hs = SVector{3,Float64}(Hb),#system momentum
    Te = SVector{3,Float64}(zeros(3)),
    Ti = SVector{3,Float64}(zeros(3)),
    eci_to_ecef = SMatrix{3,3,Float64}(eci_to_ecef),
    r_eci = SVector{3,Float64}(r0),
    r_ecef = SVector{3,Float64}(r_ecef),
    lla = SVector{3,Float64}(lla),
    v_eci = SVector{3,Float64}(v0),
    v_ecef = SVector{3,Float64}(eci_to_ecef*v0),
    rmag = Float64(rmag0)
)
    
x_sun = ComponentArray(
    num_css_detect = p.fsw.ac.ad.sun.css_detect.many_css,
    sun_vector_desired = SVector{3,Float64}(zeros(3)),
    sun_vector_measured = SVector{3,Float64}(zeros(3)),
)
x_ad = ComponentArray(
    sun = x_sun,    
    q_i2b = SVector{4,Float64}([-0.021179442940975148, 0.5496477692711645, 0.16029818602936372, 0.8195994463685526]),
    ω_i2b = SVector{3,Float64}(0,-0.0011,0),
)

x_ctrl = ComponentArray(
    attitude_error = SVector{3,Float64}(zeros(3)),
    attitude_error_lim = SVector{3,Float64}(zeros(3)), 
    attitude_error_mag = Float64(0),   
    integral_error = SVector{3,Float64}(zeros(3)),
    integral_error_lim = SVector{3,Float64}(zeros(3)),
    rate_error = SVector{3,Float64}(zeros(3)),
    rate_error_lim = SVector{3,Float64}(zeros(3)),
    torque_command = SVector{3,Float64}(zeros(3)),
    torque_command_b = SVector{3,Float64}(zeros(3)),    
    torque_gyroscopic_ffwd = SVector{3,Float64}(zeros(3)),
    torque_oci_ffwd = SVector{3,Float64}(zeros(3)),
)
x_dyn = (
    Hs_b = SVector{3,Float64}(zeros(3)),
    Hs_b_mag = Float64(0),
    Hs_b_prev = SVector{3,Float64}(zeros(3)),
    ΔHs_b = SVector{3,Float64}(zeros(3)),
    ΔHs_b_mag = Float64(0),
)

x_logic = ComponentArray(
    gnc_mode = Int64(p.fsw.ac.logic.gnc_mode.MissionScience),
    gnc_mode_prev = Int64(p.fsw.ac.logic.gnc_mode.MissionScience),
)

q_nadir = SVector{4,Float64}([-0.021179442940975148, 0.5496477692711645, 0.16029818602936372, 0.8195994463685526])
x_tgt = ComponentArray(
    r_sc_to_site_eci = SVector{3,Float64}(614519.7575155823,145081.6883704655, 271380.5845309013),
    R_eci_to_brf = SMatrix{3,3,Float64}([0.34432159624611963 0.23957867587376439 -0.9077690765966298; -0.2859308764863784 0.9477203401333472 0.14166753604571566; 0.8941549522412857 0.2111006335425908 0.39487142704974343]),
    q_ref = q_nadir,
    q_nadir =  q_nadir,
    q_yaw_steering = SA[0,0,0,1],
    ω_ref = ω,
    torque_oci_ffwd = SVector{3,Float64}(zeros(3))
)

x_ac = ComponentArray(
    ad = x_ad,
    ctrl = x_ctrl,
    dyn = x_dyn,
    logic = x_logic,
    tgt = x_tgt,
)

x_output_rw = ComponentArray(
    torque_command = SVector{4,Float64}(zeros(4)),
    torque_command_lim = SVector{4,Float64}(zeros(4)),
    current_command = SVector{4,Float64}(zeros(4)),
    rwMomRdTq = SVector{4,Float64}(zeros(4)),
    momAdj = Float64(0),
    prevMomAdj = Float64(0),
    posMomAdj =  Float64(0),
    negMomAdj = Float64(0),
    shiftDirCnt = Int64(0),
)
x_output = ComponentArray(
    rw = x_output_rw,
)
x_fsw = ComponentArray(
    ac = x_ac,
    output = x_output
)

wheel_speeds = 100*ones(4)#zeros(4) #rad/sec
x_rw = ComponentArray(
    ω = SVector{4,Float64}(zeros(4)), #wheel speed
    Tw = SVector{4,Float64}(zeros(4)), #wheel torque
    Tb = SMatrix{3,4,Float64}(zeros(3,4)), #wheel torque in body frame
    Hw = SVector{4,Float64}(zeros(4)), # wheel momentum
    Hb = SMatrix{3,4,Float64}(zeros(3,4)) # wheel momentum in body frame
)

x_mtb = ComponentArray(
    I = SVector{3,Float64}(zeros(3)),
    Mm = SVector{3,Float64}(zeros(3)),
    Mb = SVector{3,Float64}(zeros(3)),
    Tb = SVector{3,Float64}(zeros(3))
)

x_actuators = ComponentArray(
    rw = x_rw,
    mtb = x_mtb
)

x_gravity = ComponentArray(
    a = g,
    gradient_torque_b = SVector{3,Float64}(zeros(3)),
    )
x_geomagnetism = ComponentArray(
    B_ecef = SVector{3,Float64}(zeros(3)),
    B_eci = SVector{3,Float64}(zeros(3)),
    B_b = SVector{3,Float64}(zeros(3)))
x_atmosphere = ComponentArray(ρ = Float64(0.0))

x_environments = ComponentArray(
    geomagnetism = x_geomagnetism,
    gravity = x_gravity,
    atmosphere = x_atmosphere
)

x0 = ComponentArray(
    actuators = x_actuators,    
    body = x_body,
    environments = x_environments,
    fsw = x_fsw,    
    orbit = x_orbit,    
)

#used for testing and benchmarking callback functions in REPL
S = ComponentArray(
    u = x0,
    p = initModelParams(p)
)
""" ODE """
#sol = simulate(x0,p,(0,1000))