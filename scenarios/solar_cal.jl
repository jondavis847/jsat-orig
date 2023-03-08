includet("..\\src\\jsat.jl")
includet("..\\pace_parameters.jl")

p = pace_parameters()
x0 = InitialConditions(
    Mode = NominalValue(p.fsw.ac.logic.gnc_mode.RateNull),
    t = NominalValue(datetime2julian(DateTime(2024,1,9,6,43,12))),
    r_eci = NominalValue([-3.929273771183811e6,5.71264018383781e6,1.3119944637514637e6]),
    v_eci = NominalValue([84.55513412044282,1749.4937754157777,-7311.912202615174]),
    q_i2b = NominalValue([-0.021179442940975148, 0.5496477692711645, 0.16029818602936372, 0.8195994463685526]), ## CHANGE QUAT
    ω_b = DispersedValue([0, 0.5*pi/180, 0], Normal.(zeros(3),[1,2,2]*pi/180/3)),
    rwa_ω = NominalValue(zeros(4))
)

#sol = simulate(x0,p,(0,1000))