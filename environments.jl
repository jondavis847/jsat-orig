using UnPack, SatelliteToolbox
includet("orbit.jl")

"""Environments """

function environments_cb!(integrator)
   #gravity_cb!(integrator)
   #geomagnetism_cb!(integrator)
   stb_gravity_cb!(integrator)
   stb_geomagnetism_cb!(integrator)
end

function stb_gravity_cb!(S)
    S.u.environments.gravity.a = S.u.body.eci_to_ecef' * compute_g(S.p.environments.gravity.egm96_coefs, S.u.body.r_ecef)        
end
#over load length of gravity model to support component arrays
function Base.length(in::GravityModel_Coefs{Float64})
    return 0
end

function gravity_cb!(S)
    @unpack μ, R, J2, J3, J4, J5, J6 = S.p.environments.gravity
    x, y, z = S.u.body.r
    r = S.u.body.rmag

    aj2 = SVector{3}(-(3 / 2) * J2 * (μ / r^2) * (R / r)^2 * [
                         (1 - 5 * (z / r)^2) * x / r
                         (1 - 5 * (z / r)^2) * y / r
                         (3 - 5 * (z / r)^2) * z / r
                     ])

    aj3 = SVector{3}(-(1 / 2) * J3 * (μ / r^2) * (R / r)^3 * [
                         5 * (7 * (z / r)^3 - 3 * (z / r)) * x / r
                         5 * (7 * (z / r)^3 - 3 * (z / r)) * y / r
                         3 * (10 * (z / r)^2 - 35 / 3 * (z / r)^4 - 1)
                     ])

    aj4 = SVector{3}(-(5 / 8) * J4 * (μ / r^2) * (R / r)^4 * [
                         (3 - 42 * (z / r)^2 + 63 * (z / r)^4) * x / r
                         (3 - 42 * (z / r)^2 + 63 * (z / r)^4) * y / r
                         -(15 - 70 * (z / r)^2 + 63 * (z / r)^4) * z / r
                     ])

    aj5 = SVector{3}(-(1 / 8) * J5 * (μ / r^2) * (R / r)^5 * [
                         3 * (35 * (z / r) - 210 * (z / r)^3 + 231 * (z / r)^5) * x / r
                         3 * (35 * (z / r) - 210 * (z / r)^3 + 231 * (z / r)^5) * y / r
                         (15 - 315 * (z / r)^2 + 945 * (z / r)^4 - 693 * (z / r)^6)
                     ])

    aj6 = SVector{3}((1 / 16) * J6 * (μ / r^2) * (R / r)^6 * [
                         (35 - 945 * (z / r)^2 + 3465 * (z / r)^4 - 3003 * (z / r)^6) * x / r
                         (35 - 945 * (z / r)^2 + 3465 * (z / r)^4 - 3003 * (z / r)^6) * y / r
                         (2205 * (z / r)^2 - 4851 * (z / r)^4 + 3003 * (z / r)^6 - 315) * z / r
                     ])

    S.u.environments.gravity.a = aj2 + aj3 + aj4 + aj5 + aj6
end


""" Geomagnetism """
function stb_geomagnetism_cb!(S)    
    S.u.environments.geomagnetism.B = igrf(2023.0, S.u.body.lla[3], S.u.body.lla[1], S.u.body.lla[2], Val(:geodetic))
end

#need to overload EOP length to work with component arrays
function Base.length(in::EOPData_IAU1980)
    return 0
end

function geomagnetism_cb!(integrator)
    # IGRFMAG Uses International Geomagnetic Reference Field
    #
    #  Calculates the Earth magnetic field and secular variation at a specific
    #  location and time using different generations of International Geomagnetic
    #  Reference Field. 
    #  
    #  Same algorithm as that in OM, but this uses a full model (13th degree) 
    #  of the IGRF-12 coefficients, which are embedded here, not an upload table.
    #
    #  Inputs required by IGRF are:
    #   LLA    : [LAT, LON, ALT]
    #   (LAT)  : Scalar geodetic latitude in radian where north latitude is
    #          positive and south latitude is negative.
    #   (LON)  : Scalar geodetic longitude in radian where east longitude
    #          is positive and  west longitude is negative.
    #   (ALT)  : Scalar altitude in meters.
    #   DECYR  : Scalar decimal year.  Decimal year is the desired year in
    #          a decimal format that includes any fraction of the year that has
    #          already passed.
    #
    #  Output calculated for the Earth magnetic field and secular variation include:
    #   bVecECEF :Magnetic field vector in nanotesla (nT) in ECEF frame.
    #   bVecNED  :Magnetic field vector in nanotesla (nT) in NED frame.
    # 
    #   Limitations:
    #
    #   This function is valid between the heights of -1000 meters to 600000
    #   meters. 
    #   Extrapolation for out-of-range height.
    #
    #   The 13th generation is valid between the years of 1900 and 2025.
    #   Extrapolation for out-of-range height.
    # 
    #   This function has the limitations of the International Geomagnetic
    #   Reference Field (IGRF). For more information, see the IGRF web site,
    #   http://www.ngdc.noaa.gov/IAGA/vmod/igrfhw.html.
    #

    #   Modification based on MATLAB\R2018b\toolbox\aero\aero\igrfmagm.m.

    #   Reference:
    #   [1] The IGRF-13 can be found on the web at
    #       http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    #   [2] Blakely, R. J., "Potential Theory in Gravity & Magnetic Applications",
    #       Cambridge University Press, Cambridge UK, 1996.
    
    lla = ecef_to_lla(eci_to_ecef(integrator.u.orbit.epoch)*integrator.u.body.r) 
    latDeg = lla[1] * 180 / pi # [rad]->[deg]
    lonDeg = lla[2] * 180 / pi # [rad]->]deg]
    altM = lla[3] # [m]
    # IGRF max degree=13; 
    nMax = 13

    ## Load Schmidt semi-normalised spherical harmonic coefficients
    # Those are obtained from upload table.
    # igrfdata = webread('https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt');
    minYear = 2020
    gh = [-29404.8, -1450.9, 4652.5, -2499.6, 2982.0, -2991.6, 1677.0, -734.6, 1363.2, -2381.2, -82.1, 1236.2, 241.9, 525.7, -543.4, 903.0, 809.5, 281.9, 86.3, -158.4, -309.4, 199.7, 48.0, -349.7, -234.3, 363.2, 47.7, 187.8, 208.3, -140.7, -121.2, -151.2, 32.3, 13.5, 98.9, 66.0, 65.5, -19.1, 72.9, 25.1, -121.5, 52.8, -36.2, -64.5, 13.5, 8.9, -64.7, 68.1, 80.6, -76.7, -51.5, -8.2, -16.9, 56.5, 2.2, 15.8, 23.5, 6.4, -2.2, -7.2, -27.2, 9.8, -1.8, 23.7, 9.7, 8.4, -17.6, -15.3, -0.5, 12.8, -21.1, -11.7, 15.3, 14.9, 13.7, 3.6, -16.5, -6.9, -0.3, 2.8, 5.0, 8.4, -23.4, 2.9, 11.0, -1.5, 9.8, -1.1, -5.1, -13.2, -6.3, 1.1, 7.8, 8.8, 0.4, -9.3, -1.4, -11.9, 9.6, -1.9, -6.2, 3.4, -0.1, -0.2, 1.7, 3.6, -0.9, 4.8, 0.7, -8.6, -0.9, -0.1, 1.9, -4.3, 1.4, -3.4, -2.4, -0.1, -3.8, -8.8, 3.0, -1.4, 0.0, -2.5, 2.5, 2.3, -0.6, -0.9, -0.4, 0.3, 0.6, -0.7, -0.2, -0.1, -1.7, 1.4, -1.6, -0.6, -3.0, 0.2, -2.0, 3.1, -2.6, -2.0, -0.1, -1.2, 0.5, 0.5, 1.3, 1.4, -1.2, -1.8, 0.7, 0.1, 0.3, 0.8, 0.5, -0.2, -0.3, 0.6, -0.5, 0.2, 0.1, -0.9, -1.1, 0.0, -0.3, 0.5, 0.1, -0.9, -0.9, 0.5, 0.6, 0.7, 1.4, -0.3, -0.4, 0.8, -1.3, 0.0, -0.1, 0.8, 0.3, 0.0, -0.1, 0.4, 0.5, 0.1, 0.5, 0.5, -0.4, -0.5, -0.4, -0.4, -0.6]
    sv = [5.7, 7.4, -25.9, -11.0, -7.0, -30.2, -2.1, -22.4, 2.2, -5.9, 6.0, 3.1, -1.1, -12.0, 0.5, -1.2, -1.6, -0.1, -5.9, 6.5, 5.2, 3.6, -5.1, -5.0, -0.3, 0.5, 0.0, -0.6, 2.5, 0.2, -0.6, 1.3, 3.0, 0.9, 0.3, -0.5, -0.3, 0.0, 0.4, -1.6, 1.3, -1.3, -1.4, 0.8, 0.0, 0.0, 0.9, 1.0, -0.1, -0.2, 0.6, 0.0, 0.6, 0.7, -0.8, 0.1, -0.2, -0.5, -1.1, -0.8, 0.1, 0.8, 0.3, 0.0, 0.1, -0.2, -0.1, 0.6, 0.4, -0.2, -0.1, 0.5, 0.4, -0.3, 0.3, -0.4, -0.1, 0.5, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ## Interpolate/Extrapolate coefficients
    # Interpolate/Extrapolate coefficients if necessary
    # gha = zeros(ghArraySize,1); ##ok<PREALL>
    # Extrapolation
    deltaTime = decyear(integrator.u.orbit.epoch) - minYear
    # Adjusted Schmidt semi-normalised spherical harmonic coefficients
    gha = gh + deltaTime * sv

    # Calculates field components from spherical harmonic
    # WGS84
    Re = 6378.1370 # [km]
    a2 = 4.068063159076899e7
    b2 = 4.040829998408706e7

    # preallocate vectors for calculations
    sl = zeros(14)
    cl = zeros(14)
    p = zeros(119)
    q = zeros(119)

    # initialize outputs
    x = 0
    y = 0
    z = 0

    # initialize counters
    l = 1
    n = 0
    m = 1

    # convert to kilometers from meters
    elevKM = altM * 0.001

    slat_gd = sind(latDeg)

    # near polar region
    if ((90.0 - latDeg) < 0.001)
        #  300 ft. from North pole
        latDegCrt = 89.999
    elseif ((90.0 + latDeg) < 0.001)
        #  300 ft. from South pole
        latDegCrt = -89.999
    else
        latDegCrt = latDeg
    end#if
    clat_gd = cosd(latDegCrt)

    rlon = lonDeg * pi / 180 # [deg]->[rad], i.e. lla(2)
    sl[1] = sin(rlon)
    cl[1] = cos(rlon)

    # convert from geodetic to geocentric coordinates
    aa = a2 * clat_gd * clat_gd
    bb = b2 * slat_gd * slat_gd
    cc = aa + bb
    dd = sqrt(cc)
    r = sqrt(elevKM * (elevKM + 2.0 * dd) + (a2 * aa + b2 * bb) / cc)
    cd = (elevKM + dd) / r
    sd = (a2 - b2) / dd * slat_gd * clat_gd / r
    slat = slat_gd * cd - clat_gd * sd
    clat = clat_gd * cd + slat_gd * sd

    ratio = Re / r
    sqrt3 = sqrt(3.0)

    p[1] = 2.0 * slat
    p[2] = 2.0 * clat
    p[3] = 4.5 * slat * slat - 1.5
    p[4] = 3.0 * sqrt3 * clat * slat
    q[1] = -clat
    q[2] = slat
    q[3] = -3.0 * clat * slat
    q[4] = sqrt3 * (slat * slat - clat * clat)

    fn = 0
    rr = 1
    npq = (nMax * (nMax + 3)) / 2
    for k = 1:Int(npq)
        if (n < m)
            m = 0
            n = n + 1
            power = n + 2
            rr = ratio^power
            fn = n
        end
        fm = m
        if (k >= 5)
            if (m == n)
                aa = sqrt((1.0 - 0.5 / fm))
                j = k - n - 1
                p[k] = (1.0 + 1.0 / fm) * aa * clat * p[j]
                q[k] = aa * (clat * q[j] + slat / fm * p[j])
                sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1]
                cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1]
            else
                aa = sqrt(fn * fn - fm * fm)
                bb = sqrt(((fn - 1.0) * (fn - 1.0)) - (fm * fm)) / aa
                cc = (2.0 * fn - 1.0) / aa
                ii = k - n
                j = k - 2 * n + 1
                p[k] = (fn + 1.0) * (cc * slat / fn * p[ii] - bb / (fn - 1.0) * p[j])
                q[k] = cc * (slat * q[ii] - clat / fn * p[ii]) - bb * q[j]
            end
        end#if
        aa = rr * gha[l]

        if (m == 0)
            x = x + aa * q[k]
            z = z - aa * p[k]
            l = l + 1
        else
            bb = rr * gha[l+1]
            cc = aa * cl[m] + bb * sl[m]
            x = x + cc * q[k]
            z = z - cc * p[k]

            if (clat > 0)
                y = y + (aa * sl[m] - bb * cl[m]) * fm * p[k] / ((fn + 1.0) * clat)
            else
                y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat
            end
            l = l + 2
        end#if
        m = m + 1
    end#for

    x_old = x
    x = x * cd + z * sd
    z = z * cd - x_old * sd

    # Vectorize Northward, Eastward and Downward components
    bVecNED = [x, y, z]

    ## Convert magnetic field from NED to ECEF
    C_ecef2ned = [
        -sind(latDeg) 0 cosd(latDeg)
        0 1 0
        -cosd(latDeg) 0 -sind(latDeg)
    ] *
                 [
        cosd(lonDeg) sind(lonDeg) 0
        -sind(lonDeg) cosd(lonDeg) 0
        0 0 1
    ]


    # Rotate to Earth-Centered, Earth-Fixed frame
    bVecECEF = C_ecef2ned' * bVecNED

    integrator.u.environments.geomagnetism.B = bVecECEF
end# igrfmag()

