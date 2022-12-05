using StaticArrays,SatelliteToolbox,Dates

# Julian date calculation

jd(Y, M, D, h, m, s) = 1721013.5 + 367 * Y - trunc(7 / 4 * (Y + trunc((M + 9) / 12))) + trunc(275 * M / 9) + D + (60 * h + m + s / 60) / 1440
jd(t::DateTime) = jd(Dates.year(t),Dates.month(t),Dates.day(t),Dates.hour(t),Dates.minute(t),Dates.second(t)+Dates.millisecond(t))
date_to_jd(t::DateTime) = SatelliteToolbox.date_to_jd(Dates.year(t),Dates.month(t),Dates.day(t),Dates.hour(t),Dates.minute(t),Dates.second(t)+Dates.millisecond(t))
#=
Base.length(in::DateTime)= 0
Base.UNITLESS_ABS2(in::DateTime) = 0

function eci_to_ecef(epoch)
    @unpack Y, M, D, h, m, s = epoch
    # Number of julian centuries elapsed from the epoch j2000.0 to zero hours of the date in question
    T0 = (jd(Y, M, D, 0, 0, 0) - 2451545) / 36525

    # Greenwich Mean Sidereal Time angle
    θGMST = ((24110.54841 + 8640184.812866 * T0 + 0.093104 * T0^2 - 6.2e-6 * T0^3 + 1.002737909350795 * (3600 * h + 60 * m + s)) % 86400) / 240 #2.70

    R_EciToEcef = @SMatrix [
        cos(θGMST) sin(θGMST) 0
        -sin(θGMST) cos(θGMST) 0
        0 0 1
    ]
    return R_EciToEcef
end
=#
# ECEF to Latitude,Longitude,Altitude (LLA)
#Markley/Crassidis equations 2.77(a-g)
function ecef_to_lla_not_working(r_ecef)
    a = 6378137.0 #m
    b = 6356752.3142 #m
    x, y, z = r_ecef
    e² = 1 - b^2 / a^2
    ϵ² = a^2 / b^2 - 1
    ρ = √(x^2 + y^2)
    p = abs(z) / ϵ²
    s = ρ^2 / (e² * ϵ²)
    q = p^2 - b^2 + s
    u = p / √(q)
    v = b^2 * u^2 / q
    P = 27 * v * s / q
    Q = (√(P + 1) + √(P))^(2 / 3)
    t = (1 + Q + 1 / Q) / 6
    c = √(u^2 - 1 + 2 * t)
    w = (c - u) / 2
    d = sign(z) * √(q) * (w + (√(t^2 + v) - u * w - t / 2 - 1 / 4)^1 / 2)
    N = a * √(1 + ϵ² * d^2 / b^2)
    λ = asin((ϵ² + 1) * (d / N)) #latitude
    h = ρ * cos(λ) + z * sin(λ) - a^2 / N #altitude
    ϕ = atan(y, x) #longitude

    return @SVector [λ, ϕ, h]
end

function ecef_to_lla(r_ecef)
    # Journal of Geodesy (2002) 76: 451–454
    # DOI:10.1007/s00190-002-0273-6
    #-----------------------------------------------------
    x, y, z = r_ecef

    # geodetic (ellipsoid) parameters (a,f) or (a,e)
    # semimajor axis, determined from a combination of Doppler satellite and astro-geodetic data 
    # flattening, determined from satellite data
    a = 6378137.0 # semimajor axis (6378137 m, WGS84 value)
    ecc = 0.0818
    b = a * sqrt(1 - ecc^2) # semiminor axis
    esq = ecc^2 # eccentricity square; ecc^2 = 2*f - f^2;
    # f = 1-b/a; # flattening (1/298.25, WGS84 value)

    # Start 15-step procedure
    r = sqrt(x^2 + y^2)
    r2 = r^2
    b2 = b^2
    a2 = a^2
    Z2 = z^2
    esq2 = esq^2
    epsq = (a2 / b2) - 1 # e_prime square
    F = 54 * b2 * Z2
    G = r2 + (1 - esq) * Z2 - esq * (a2 - b2)
    C = (esq2 * F * r2) / (G * G * G)
    S = (1 + C + sqrt(C^2 + 2 * C))^(1 / 3)
    P = F / (3 * (S + (1 / S) + 1)^2 * G^2)
    Q = sqrt(1 + 2 * esq2 * P)
    r0Det = 0.5 * a2 * (1 + 1 / Q) - (1 - esq) * P * Z2 / (Q * (1 + Q)) - 0.5 * P * r2
    r0 = -esq * P * r / (1 + Q) + sqrt(abs(r0Det))
    U = sqrt((r - esq * r0)^2 + Z2)
    V = sqrt((r - esq * r0)^2 + (1 - esq) * Z2)
    Z0 = b2 * z / (a * V)
    # Solution ...
    lat = atan((z + epsq * Z0), r)
    lon = atan(y, x)
    alt = U * (1 - b2 / (a * V))
    return @SVector [lat, lon, alt]
end

function decyear(datetime)    
    leapyear = Dates.isleapyear(datetime)
    ndays = leapyear ? 366 : 365
    dayofyear = Dates.Day(floor(datetime,Dates.Day) - floor(datetime,Dates.Year))    
    decyear = Dates.year(datetime) + dayofyear.value/ndays
    return decyear
end