using Interpolations, Distributions, SatelliteToolbox

struct ModelParameter
    name
    value
    dist
end
Base.length(in::ModelParameter) = length(in.value)

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
    b_to_rw=SMatrix{4,3,Float64}(b_to_rw),
    rw_to_b=SMatrix{3,4,Float64}(b_to_rw'),
    EnableMinimax=Int64(0),       # Flag to select distribution logic (0=DIS, 1=ENA)
    momRdGain=Float64(0.0),              # Momentum redistribution gain
    nullVec= SVector{4,Float64}([1, -1, 1, -1]),     # Null vector of mounting matrix, RwaToBcs
    whRdTqLim=Float64(0.225),            # Wheel torque limit for torque redistribution [Nm]
    cmdNullMomZeroAvoid=Float64(0),        # Null momentum to avoid particular wheel speeds that cause resonance [Nm]
    rwaZeroMomThrd=Float64(0),            # Wheel speed is considered close to zero if wheel momentum magnitude is below this value [Nms]
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
    EnableMomMax=0   # Flag to select enforcement of wheel momentum magnitude limit (0=DIS, 1=ENA)
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

p_fsw = (
    ac = p_ac,
)

wheel_km = 0.077 * ones(4)
wheel_J = 0.231 * ones(4)
wheel_trq_max = 0.40 * ones(4)
wheel_momentum_max = 70 * ones(4)

p_rw = (
    km=SVector{4,Float64}(wheel_km), #motor constant
    J=SVector{4,Float64}(wheel_J), #wheel inertia
    a=SMatrix{3,4,Float64}([a'...;]'), #axis of rotation in reference frame, convert to matrix since component arrays can't have arrays of arrays :(  
    WhTqMax=wheel_trq_max,   # Max wheel torque (N*m)
    WhMomMax=wheel_momentum_max,     # Max wheel momentum (N*m*s)    
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
