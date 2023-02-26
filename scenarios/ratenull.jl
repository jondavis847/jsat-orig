includet("..\\src\\jsat.jl")
includet("..\\pace_parameters.jl")

p = pace_parameters()

x0 = InitialConditions(
    Mode = NominalValue(p.fsw.ac.logic.gnc_mode.RateNull),
    t = NominalValue(datetime2julian(DateTime(2023,1,1,12,0,0))),
    r_eci = NominalValue([-6.32053381835774e6, -1.49222035669749e6, -2.77429961381375e6]),
    v_eci = NominalValue([2.56786681792351e3, 1.79527532306182e3, -6.82715713742553e3]),
    q_i2b = NominalValue([-0.021179442940975148, 0.5496477692711645, 0.16029818602936372, 0.8195994463685526]),
    ω_b = DispersedValue([0, 0.5*pi/180, 0], Normal.(zeros(3),[1,2,2]*pi/180/3)),
    rwa_ω = NominalValue(zeros(4))
)

#sol = simulate(x0,p,(0,1000))