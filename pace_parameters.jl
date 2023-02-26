includet("src\\model.jl")
using Interpolations, Distributions, SatelliteToolbox,ComponentArrays,LinearAlgebra

function pace_parameters()
"""Config"""
p_config = ComponentArray(
    geomagnetism=true,
    atmosphere=true,    
    montecarlo = false,
)

"""Mass Properties"""

mass = DispersedValue(1506.27, Normal(1506.27, 1506.27 * 0.05 / 3))

inertia = DispersedValue(
    # Nominal Value
    [
        1529.097 -58.0883 -26.71023
        -58.0883 1400.682 83.51491
        -26.71023 83.51491 2320.778
    ], Normal.(
        # Mean Value    
        [
            1529.097 -58.0883 -26.71023
            -58.0883 1400.682 83.51491
            -26.71023 83.51491 2320.778
        ],
        # Sigma Values
        [
            1529.097*0.1 100 100
            100 1400.682*0.1 100
            100 100 2320.778*0.1
        ] ./ 3))

        
cg = DispersedValue([-1.311708,-0.122079,-0.0302493], 
    Normal.([-1.311708,-0.122079,-0.0302493],abs.([-1.311708,-0.122079,-0.0302493]*0.1/3)))

p_body = ComponentArray(
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

css_detect = ComponentArray(
    no_css=Int64(0),
    one_css=Int64(1),
    two_css=Int64(2),
    many_css=Int64(3)
)
p_sun = ComponentArray(
    css_detect=css_detect,
    sun_vector_desired=SVector{3,Float64}(0, 0, -1),
    rates_desired=SVector{3,Float64}(zeros(3))
)

p_ad = ComponentArray(
    sun=p_sun,
)

p_ctrl = ComponentArray(
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

gnc_mode = ComponentArray(
    LaunchMode=0,
    SunSafe=1,
    RateNull=2,
    MissionScience=3,
    InertialReference=4,
    LunarCal=5,
    DeltaV=6,
    DeltaH=7,
)

p_logic = ComponentArray(
    gnc_mode=gnc_mode,
)

p_tgt = ComponentArray(
    orbital_rate=SVector{3,Float64}(0, -0.0011, 0),
)
p_ac = ComponentArray(
    ad=p_ad,
    ctrl=p_ctrl,
    logic=p_logic,
    tgt=p_tgt
)

wheel_km = 0.077 * ones(4)
wheel_J = 0.231 * ones(4)
wheel_trq_max = 0.40 * ones(4)
wheel_momentum_max = 70 * ones(4)

p_output_rw = ComponentArray(
    km=SVector{4,Float64}(wheel_km), #motor constant
    J=SVector{4,Float64}(wheel_J), #wheel inertia
    WhTqMax=wheel_trq_max,   # Max wheel torque (N*m)
    WhMomMax=wheel_momentum_max,     # Max wheel momentum (N*m*s)   
    b_to_rw=SMatrix{4,3,Float64}(b_to_rw),
    rw_to_b=SMatrix{3,4,Float64}(b_to_rw'),
    EnableMinimax=Int64(0),       # Flag to select distribution logic (0=DIS, 1=ENA)
    momRdGain=Float64(0),              # Momentum redistribution gain
    nullVec=SVector{4,Float64}([1, -1, 1, -1]),     # Null vector of mounting matrix, RwaToBcs
    whRdTqLim=Float64(0.225),            # Wheel torque limit for torque redistribution [Nm]
    cmdNullMomZeroAvoid=Float64(0),       # Null momentum to avoid particular wheel speeds that cause resonance [Nm]
    rwaZeroMomThrd=Float64(0 * 20 * 2 * pi / 60),            # Wheel speed is considered close to zero if wheel momentum magnitude is below this value [Nms]
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

p_output = ComponentArray(
    rw=p_output_rw,
)

p_fsw = ComponentArray(
    ac=p_ac,
    output=p_output,
)


p_rw = ComponentArray(
    km=DispersedValue(wheel_km, Normal.(wheel_km,0.05.*wheel_km./3)), #motor constant
    J=DispersedValue(wheel_J, Normal.(wheel_J,0.05.*wheel_J./3)), #wheel inertia
    a=SMatrix{3,4,Float64}(hcat(a...)), #axis of rotation in reference frame, convert to matrix since component arrays can't have arrays of arrays :(       
    f_bemf = SVector{4,Float64}(zeros(4)),
    knee = DispersedValue(53.052*ones(4),Normal.(53.052,0.05*53.052*ones(4)/3)),
    H_m = DispersedValue(3.75e-4*ones(4), Normal.(3.75e-4*ones(4),0.05*3.75e-4*ones(4)./3)), #momentum curve slope
    H_b = DispersedValue(0.505*ones(4),Normal.(0.505*ones(4),0.05*0.505/3*ones(4))), #momentum curve bias
    T_m = DispersedValue(0.02304*ones(4),Normal.(0.02304*ones(4),0.05*0.02304*ones(4)/3)),
    T_b = DispersedValue(1.70742857142857*ones(4),Normal.(1.70742857142857*ones(4),0.05*1.70742857142857*ones(4)/3)), #torque curve bias
    n_poles = 8,
    n_phases = 3,
    cogging_amplitude = 0.008,
    cogging_phase = 0,
    ripple_coeff = 0.016,
    ripple_phase = 0,
    f_coulomb = DispersedValue(0.02.*ones(4),Normal.(0.02.*ones(4),0.02*0.05/3 .*ones(4))),
    f_stiction = DispersedValue(0.03.*ones(4),Normal.(0.03.*ones(4),0.03*0.05/3 .*ones(4))),
    f_viscous = DispersedValue(8e-5.*ones(4),Normal.(8e-5.*ones(4),8e-5*0.05/3 .*ones(4))),
    f_windage = DispersedValue(3.3e-8.*ones(4),Normal.(3.3e-8.*ones(4),3.3e-8*0.05/3 .*ones(4))),
    beta_s = 20,
    stiction_speed = DispersedValue(ones(4),Normal.(ones(4),0.5/3 .*ones(4))),#1e5#arbitrary  
    R_motor = DispersedValue(2.025.*ones(4),Normal.(2.025 .*ones(4),2.025*0.05/3).*ones(4)),
)

current = Float64[-450, -405, -360, -315, -270, -225, -180, -90, 0, 90, 180, 225, 270, 315, 360, 405, 450]
moment = Float64[-450, -440, -430, -420, -400, -355, -300, -150, 0, 150, 300, 355, 400, 420, 430, 440, 450]
p_mtb = ComponentArray(
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

p_actuators = ComponentArray(rw=p_rw, mtb=p_mtb)

egm96_model = parse_icgem("EGM96.gfc")
egm96_coefs = create_gravity_model_coefs(egm96_model)
egm96_coefs = load_gravity_model(EGM96())

p_gravity = ComponentArray(
    egm96_coefs=egm96_coefs,
)

eop_IAU1980 = get_iers_eop()
p_geomagnetism = ComponentArray(
    eop_IAU1980=eop_IAU1980,
)

p_environments = ComponentArray(
    gravity=p_gravity,
    geomagnetism=p_geomagnetism
)
p = ComponentArray(
    config=p_config,
    body=p_body,
    fsw=p_fsw,
    actuators=p_actuators,
    environments=p_environments,
)

end