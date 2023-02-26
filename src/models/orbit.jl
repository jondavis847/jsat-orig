using StaticArrays,SatelliteToolbox,Dates

# Julian date calculation

jd(Y, M, D, h, m, s) = 1721013.5 + 367 * Y - trunc(7 / 4 * (Y + trunc((M + 9) / 12))) + trunc(275 * M / 9) + D + (60 * h + m + s / 60) / 1440
jd(t::DateTime) = jd(Dates.year(t),Dates.month(t),Dates.day(t),Dates.hour(t),Dates.minute(t),Dates.second(t)+Dates.millisecond(t))
date_to_jd(t::DateTime) = SatelliteToolbox.date_to_jd(Dates.year(t),Dates.month(t),Dates.day(t),Dates.hour(t),Dates.minute(t),Dates.second(t)+Dates.millisecond(t))

function decyear(datetime)    
    leapyear = Dates.isleapyear(datetime)
    ndays = leapyear ? 366 : 365
    dayofyear = Dates.Day(floor(datetime,Dates.Day) - floor(datetime,Dates.Year))    
    decyear = Dates.year(datetime) + dayofyear.value/ndays
    return decyear
end

function orbit_cb!(S)
    time!(S)
    sunVector!(S)
    moonVector!(S)
end

function time!(S)
    S.u.orbit.epoch_prev = S.u.orbit.epoch
    S.u.orbit.epoch = S.u.orbit.epoch + (S.t - S.tprev) / 86400.0
    return nothing
end


function sunVector!(S)
    sunMOD = sun_position_i(S.u.orbit.epoch)
    mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),S.u.orbit.epoch)
    S.u.orbit.sun_position_eci = mod_to_j2000 * sunMOD
    S.u.orbit.sun_vector_eci = normalize(S.u.orbit.sun_position_eci - S.u.body.r_eci)
    S.u.orbit.sun_vector_b = qvrot(S.u.body.q,S.u.orbit.sun_vector_eci)
end

function moonVector!(S)
    moonMOD = moon_position_i(S.u.orbit.epoch)
    mod_to_j2000 = r_eci_to_eci(MOD(),J2000(),S.u.orbit.epoch)
    S.u.orbit.moon_position_eci = mod_to_j2000 * moonMOD
    S.u.orbit.moon_vector_eci = normalize(S.u.orbit.moon_position_eci - S.u.body.r_eci)
    S.u.orbit.moon_vector_b = qvrot(S.u.body.q,S.u.orbit.moon_vector_eci)
end
