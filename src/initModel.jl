includet("src\\jsat.jl")
includet("pace_parameters.jl")

using Interpolations, Distributions, SatelliteToolbox
function initModel(x::States,p::Parameters)



struct ModelParameter
    name
    value
    dist
end
Base.length(in::ModelParameter) = length(in.value)
Base.one(in::ModelParameter) = 1

"""Config"""


p_config = (
    geomagnetism=true,
    atmosphere=true,
)

"""Mass Properties"""

mass = ModelParameter("Mass", 1506.27, Normal(1506.27, 1506.27 * 0.05 / 3))
ixx = ModelParameter("Ixx", 1529.097, Normal(1529.097, 1529.097 * 0.1 / 3))
iyy = ModelParameter("Iyy", 1400.682, Normal(1400.682, 1400.682 * 0.1 / 3))
izz = ModelParameter("Izz", 2320.778, Normal(2320.778, 2320.778 * 0.1 / 3))
ixy = ModelParameter("Ixy", -58.0883, Normal(0, 100 / 3))
ixz = ModelParameter("Ixz", -26.71023, Normal(0, 100 / 3))
iyz = ModelParameter("Iyz", 83.51491, Normal(0, 100 / 3))

inertia = [
    ixx ixy ixz
    ixy iyy iyz
    ixz iyz izz
]

cg_x = ModelParameter("cg_x", -1.311708, Normal(-1.311708, 0.01 / 3))
cg_y = ModelParameter("cg_y", -0.122079, Normal(-0.122079, 0.01 / 3))
cg_z = ModelParameter("cg_z", -0.0302493, Normal(-0.0302493, 0.01 / 3))

cg = [cg_x, cg_y, cg_z]

p_body = (
    J=inertia,
    mass=mass,
    cg=cg,
)

a1 = [0.76604, 0.64279, 0]
a2 = [0, 0.64279, 0.76604]
a3 = [-0.76604, 0.64279, 0]
a4 = [0, 0.64279, -0.76604]
a = [a1, a2, a3, a4]

b_to_rw = [a'...;]

css_detect = (
    no_css=Int64(0),
    one_css=Int64(1),
    two_css=Int64(2),
    many_css=Int64(3)
)
p_sun = (
    css_detect = css_detect,
    sun_vector_desired = SVector{3,Float64}(0,0,-1),
    rates_desired = SVector{3,Float64}(zeros(3))
)

p_ad = (
    sun = p_sun,
)

p_ctrl = (
    Kd_RN=SVector{3,Float64}([0.7000, 0.5372, 0.7000]),
    Kp_SS=SVector{3,Float64}(0.001 * ones(3)),
    Kd_SS=SVector{3,Float64}(0.05 * ones(3)),
    Kp_MS=SVector{3,Float64}([0.2853, 0.0288, 0.2853]),
    Kd_MS=SVector{3,Float64}([0.7478, 0.2376, 0.7478]),
    Ki_MS=SVector{3,Float64}([0.0071, 0.0010, 0.0071]),
    aeLim_MS=SVector{3,Float64}([0.25, 0.75, 0.5] * pi / 180),
    reLim_MS=SVector{3,Float64}(1e6 * ones(3)),
    iaeLim_MS=SVector{3,Float64}(0.0333 * ones(3)),
    Kp_IR=SVector{3,Float64}([0.2853, 0.0900, 0.2853]),
    Kd_IR=SVector{3,Float64}([0.5609, 0.2700, 0.6561]),
    Ki_IR=SVector{3,Float64}([0.0071, 0.0010, 0.0071]),
    aeAngLim_nonlin_IR=Float64(5 * 180 / pi),
    aeLim_IR=SVector{3,Float64}([0.25, 0.75, 0.5] * pi / 180),
    reLim_IR=SVector{3,Float64}(1e6 * ones(3)),
    iaeLim_IR=SVector{3,Float64}(0.0333 * ones(3)),
    Kp_LC=SVector{3,Float64}([0.2500 0.4000 0.2500]),
    Ki_LC=SVector{3,Float64}([0.0000 0.0000 0.0000]),
    aeLim_LC=SVector{3,Float64}([0.25, 0.75, 0.5] * pi / 180),
    reLim_LC=SVector{3,Float64}(1e6 * ones(3)),
    Kp_DV=SVector{3,Float64}([0.2823, 0.0288, 0.2853]),
    Kd_DV=SVector{3,Float64}([0.7477, 0.2375, 0.7477]),
    Ki_DV=SVector{3,Float64}([0.0000, 0.0000, 0.0000]),
    Kp_DH=SVector{3,Float64}([0.2823, 0.0288, 0.2853]),
    Kd_DH=SVector{3,Float64}([0.7477, 0.2375, 0.7477]),
    moi=SMatrix{3,3,Float64}(diagm([1815.07, 1972.58, 2589.39])), # mean MoI [kg-m^2]
    unitQuat_tol=Float64(0), # unit quat check tolerance 
    refrate=SVector{3,Float64}([0, -0.0011, 0]),    
)

gnc_mode = (
    LaunchMode=0,
    SunSafe=1,
    RateNull=2,
    MissionScience=3,
    InertialReference=4,
    LunarCal=5,
    DeltaV=6,
    DeltaH=7,
)

p_logic = (
    gnc_mode = gnc_mode,
    )

p_tgt = (
    orbital_rate = SVector{3,Float64}(0,-0.0011,0),
)
p_ac = (
    ad=p_ad,
    ctrl=p_ctrl,
    logic=p_logic,
    tgt=p_tgt
)

wheel_km = 0.077 * ones(4)
wheel_J = 0.231 * ones(4)
wheel_trq_max = 0.40 * ones(4)
wheel_momentum_max = 70 * ones(4)

p_output_rw = (
    km=SVector{4,Float64}(wheel_km), #motor constant
    J=SVector{4,Float64}(wheel_J), #wheel inertia
    WhTqMax=wheel_trq_max,   # Max wheel torque (N*m)
    WhMomMax=wheel_momentum_max,     # Max wheel momentum (N*m*s)   
    b_to_rw=SMatrix{4,3,Float64}(b_to_rw),
    rw_to_b=SMatrix{3,4,Float64}(b_to_rw'),
    EnableMinimax=Int64(0),       # Flag to select distribution logic (0=DIS, 1=ENA)
    momRdGain=Float64(0),              # Momentum redistribution gain
    nullVec= SVector{4,Float64}([1, -1, 1, -1]),     # Null vector of mounting matrix, RwaToBcs
    whRdTqLim=Float64(0.225),            # Wheel torque limit for torque redistribution [Nm]
    cmdNullMomZeroAvoid= Float64(0),       # Null momentum to avoid particular wheel speeds that cause resonance [Nm]
    rwaZeroMomThrd=Float64(0*20*2*pi/60),            # Wheel speed is considered close to zero if wheel momentum magnitude is below this value [Nms]
    shiftDirThrd=Int64(1e3),             # Threshold for persistence counter value to exceed before actually changing momentum adjustment direction
    uijUD=[    # upper diagonal terms of uij in Linf computation
        0.5077 0.0000 -0.5077 -0.5077 -1.0154 -0.5077
        -0.8520 0.0000 0.8520 -0.8520 0.0000 -0.8520
        0.5077 1.0154 0.5077 0.5077 0.0000 -0.5077
    ],
    wijUD=[# upper diagonal terms of wij in Linf computation
        0.3264 0.0000 -0.3264 -0.3264 -0.6527 -0.3264
        -0.3889 0.0000 0.3889 -0.3889 0.0000 -0.3889
        0.3264 0.6527 0.3264 0.3264 0.0000 -0.3264
    ],
    CoulombFricCoeff=zeros(4), # Coulomb friction coefficient of each wheel [Nm]
    ViscousFricCoeff=zeros(4), # Viscous friction coefficient of each wheel [Nms/rad]
    rwaDragGain=0,          # Scalar gain to compensate wheel drag
    EnableMomMax=0,   # Flag to select enforcement of wheel momentum magnitude limit (0=DIS, 1=ENA)
)

p_output = (
    rw = p_output_rw,
)

p_fsw = (
    ac = p_ac,
    output = p_output,
)


p_rw = (
    km=SVector{4,Float64}(wheel_km), #motor constant
    J=SVector{4,Float64}(wheel_J), #wheel inertia
    a=SMatrix{3,4,Float64}([a'...;]'), #axis of rotation in reference frame, convert to matrix since component arrays can't have arrays of arrays :(       
)

current = Float64[-450, -405, -360, -315, -270, -225, -180, -90, 0, 90, 180, 225, 270, 315, 360, 405, 450]
moment = Float64[-450, -440, -430, -420, -400, -355, -300, -150, 0, 150, 300, 355, 400, 420, 430, 440, 450]
current_to_moment_curve = LinearInterpolation(current, moment)
p_mtb = (
    current_to_moment=linear_interpolation(current, moment),
    current_lookup=Float64[-450, -405, -360, -315, -270, -225, -180, -90, 0, 90, 180, 225, 270, 315, 360, 405, 450],
    moment_lookup=Float64[-450, -440, -430, -420, -400, -355, -300, -150, 0, 150, 300, 355, 400, 420, 430, 440, 450],
    residual_moment=SVector{3,Float64}(1.3 * ones(3)),
    mtb_to_brf=SMatrix{3,3,Float64}([
        -1 0 0
        0 1 0
        0 0 -1
    ]),
    dipole_gain=1e5,
    dipole_bias=zeros(3),
    dipole_maxlinear=340,
    Moment2Current_Slope=[0.5800, 0.5690, 0.5690],
    Moment2Current_Intercept=[-0.0440, 0.0090, -0.0220],
)

p_actuators = (rw=p_rw, mtb=p_mtb)

egm96_model = parse_icgem("EGM96.gfc")
egm96_coefs = create_gravity_model_coefs(egm96_model)
egm96_coefs = load_gravity_model(EGM96())

p_gravity = (
    egm96_coefs=egm96_coefs,
)

eop_IAU1980 = get_iers_eop()
p_geomagnetism = (
    eop_IAU1980=eop_IAU1980,
)

p_environments = (
    gravity=p_gravity,
    geomagnetism=p_geomagnetism
)
p = (
    config=p_config,
    body=p_body,
    fsw=p_fsw,
    actuators=p_actuators,
    environments=p_environments,
)



epoch = Dates.datetime2julian(x.t0)
eci_to_ecef = r_eci_to_ecef(J2000(), ITRF(), epoch,eop_IAU1980)
r_ecef = eci_to_ecef*x.r0
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