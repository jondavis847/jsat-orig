module lunisolarEphem 
#=
    # This block predicts lunar (Moon) and solar (Sun) ephemeris 
    # (position and velocity) in ECI J2000 frame. Units are m, m/s. 
    #
    # A byproduct is True-of-date to inertial J2000 frame conversion 
    # matrix (TodToJ2000Mtx). This will be used in ECI2ECEF conversion.
    #
    # Author: Huaizu You (GSFC-5910)
    # Date: 2019.05.18 
    
    # Public, tunable properties
    properties
    =#
        SunEpoch	     = 2451545.0	            
        CosSunIncl	     = 0.91748408310236	    
        SinSunIncl	     = 0.39777249434045	    
        SunMeanCoef1	 = 6.24004076807029	    
        SunMeanCoef2	 = 1.720197e-2	        
        SolLongCoef1	 = 4.89496787343582	    
        SolLongCoef2	 = 1.720279e-2	        
        SunAnomCrCf1	 = 0.03342305517569	    
        SunAnomCrCf2	 = 3.490658503988659e-4	
        ErthSunDst0	     = 1.49597870e8	        
        ErthSunDst1	     = 1.00014	            
        ErthSunDst2	     = 0.01671	            
        ErthSunDst3	     = 0.00014	            
         
        Inclination	     = 8.98041063e-2         
        AscNodeCf0	     = 0.963504179           
        AscNodeCf1	     = -9.24217638e-4        
        LongCoef0		 = 3.18555872e-1         
        LongCoef1		 = 2.29971502e-1         
        PerigeeCf0	     = 3.367047043           
        PerigeeCf1	     = 1.94435979e-3         
        ElongCf0	     = 1.72160213            
        ElongCf1	     = 2.12768711e-1         
        SolAnomCf0	     = 6.229620611           
        SolAnomCf1	     = 1.720197e-2           
        OblEclCf0	     = 4.0909294365474e-1    
        OblEclCf1	     = 0                     
        ErthRad		     = 6.3781363e3           
        JulTimeCoef    	 = 2446065.5             
        JulTime2000	     = 2451545.0             
        J2000A1	         = 2.43817435e-2         
        J2000A2	         = 5.38608607e-6         
        J2000B1	         = 2.227870187e-4        
        J2000B2	         = -1.60570291e-7        
        J2000Cprim1    	 = 8.94240386e-2         
        J2000Cprim2    	 = -2.01648011e-2        
        J2000Cprim3	     = -3.42782665e-6 
        
        GmstCf0          = 6.697374558   
        GmstCfD0         = 0.06570982441908
        GmstCfH          = 1.00273790935
        GmstCfT2         = 0.000026
    
   
        SECS_PER_DAY     = 86400                
        DAYS_PER_CENTURY = 36525                
        DEG_TO_RAD       = pi/180               
        ARCSEC_TO_RAD    = 1/3600*pi/180        

        function setupImpl(this)
            resetImpl(this);
        end# setupImpl()

        function stepImpl(this, JD)
            sunPosTod,  sunVelTod  = solarModel(JD); # Solar ephemeris
            moonPosTod, moonVelTod = lunarModel(JD); # Lunar ephemeris
            R_TodToJ2000  = TodToJ2000Calculation(JD); # True-Of-Date to inertial J2000 direction cosine matrix
            sunPosJ2000   = R_TodToJ2000 * sunPosTod * 1e3; # [m]
            sunVelJ2000   = R_TodToJ2000 * sunVelTod * 1e3; # [m/s]
            moonPosJ2000  = R_TodToJ2000 * moonPosTod* 1e3; # [m]
            moonVelJ2000  = R_TodToJ2000 * moonVelTod* 1e3; # [m/s]
            R_EcefToEci, GHA = EcefToEciCalculation(JD, R_TodToJ2000); # Greenwich Hour Angle in [rad]
            return sunPosJ2000, moonPosJ2000, sunVelJ2000, moonVelJ2000, R_EcefToEci, GHA
        end# stepImpl()
            
        function resetImpl(this)
            this.isInit  = [true true];
            this.prevSunPosGci  = [0 0 0]';
            this.prevMoonPosGci = [0 0 0]';
            this.prevSunModelTimeSecs  = 0;
            this.prevMoonModelTimeSecs = 0;
        end# resetImpl()
        
       
        function  [sunPosGci, sunVelGci, isValid] = solarModel(this, JD)  
            ## solarModel calculates apparent Sun pos/vel in Geocentric Inertial Coordinate (True-Of-Date)
            # Input: 
            #    JD = Julian Date in days
            # Outputs:
            #    sunPosGci = Sun pos vector from Earth in TOD (km)
            #    sunVelGci = Sun vel vector in TOD (km/s)
            #    isValid   = true if (JD-JD_prev)>0
            
            # Compute time since J2000. Onboard time is calculated in seconds and referenced to 
            # May 24, 1968, 0 hour UTC and is corrected by the sun epoch to calculate J2000 time in days.
            sunModelTimeDays = JD - this.SunEpoch;
            sunModelTimeSecs = sunModelTimeDays * this.SECS_PER_DAY;
            
            # Compute the mean anomaly of the Sun in radians from the linear polynomial fit.  Put the value in the range 0 to 2*pi.
            solarMeanAnomaly = wrapInTwoPi(this.SunMeanCoef1 + this.SunMeanCoef2 * sunModelTimeDays);
            
            # Compute the mean longitude of the Sun, corrected for abberation, in radians from the linear polynomial fit.  Put the value in the range 0 to 2*pi.
            solarMeanLongitude = wrapInTwoPi(this.SolLongCoef1 + this.SolLongCoef2 * sunModelTimeDays);
            
            # Compute the ecliptic longitude by correcting the mean solar longitude.This is approximately the ‘right ascension of the sun’ with respect to Geocentric Inertial (GCI) coordinates.
            correctionEM = this.SunAnomCrCf1 * sin(  solarMeanAnomaly) ...
                         + this.SunAnomCrCf2 * sin(2*solarMeanAnomaly);
            eclipticLongitude = solarMeanLongitude + correctionEM;
            
            # Compute the unit vector from the center of the Earth to the sun's position specified in GCI..
            sunUnitVec    = [0 0 0]';
            sunUnitVec(1) = cos(eclipticLongitude);
            sunUnitVec(2) = sin(eclipticLongitude) * this.CosSunIncl;
            sunUnitVec(3) = sin(eclipticLongitude) * this.SinSunIncl;
            
            # Compute the current distance from the Earth to the Sun, and the matching vector. 
            # Buffer the previous cycle's value before calculating this cycle's sunPosGci.
            distFromEarthToSun = this.ErthSunDst0 * ...
                               ( this.ErthSunDst1 ...
                               - this.ErthSunDst2 * cos(  solarMeanAnomaly) ...
                               - this.ErthSunDst3 * cos(2*solarMeanAnomaly) );
            sunPosGci = sunUnitVec * distFromEarthToSun; # [km]
            
            # Compute the solar velocity. Requires two cycles of position and time.
            isValid = false;
            if ~this.isInit(1)
                solarDeltaTime = sunModelTimeSecs - this.prevSunModelTimeSecs; # [sec]
                sunVelGci = (sunPosGci - this.prevSunPosGci) / solarDeltaTime; # [km/s]
                if ( abs(solarDeltaTime) > 0)
                    isValid = true;
                end
            else
                sunVelGci = [0 0 0]';
                this.isInit(1) = false;
            end#if  
            
            # update states
            this.prevSunPosGci = sunPosGci;
            this.prevSunModelTimeSecs = sunModelTimeSecs;
        end# solarModel()
        
        
        function [moonPosGci, moonVelGci, isValid] = lunarModel(this,JD)
            ## lunarModel calculates apparent Moon pos/vel in Geocentric Inertial Coordinate (True-Of-Date)
            # Input: 
            #    JD = Julian Date in days
            # Outputs:
            #    moonPosGci = Moon pos vector from Earth in TOD (km)
            #    moonVelGci = Moon vel vector in TOD (km/s)
            #    isValid    = true if (JD-JD_prev)>0
            
            # Calculate time in days since the coefficients epoch.
            lunarModelTimeDays = JD - this.JulTimeCoef;
            moonModelTimeSecs = lunarModelTimeDays * this.SECS_PER_DAY;
            
            # Compute the mean longitude of the moon, measured in the ecliptic from the mean equinox of date to the mean ascending node, then along mean orbit, at the current time.
            lunarMeanLongitude = this.LongCoef0 + this.LongCoef1*lunarModelTimeDays;
            
            # Compute the mean longitude of lunar perigee measured in the ecliptic from the mean equinox of date to the mean ascending node, then along mean orbit, at the current time..
            lunarMeanPerigee = this.PerigeeCf0 + this.PerigeeCf1*lunarModelTimeDays;
            
            # Compute the longitude of the moon’s mean ascending node, measured in the ecliptic from the mean equinox of date.
            lunarMeanAscNode = this.AscNodeCf0 + this.AscNodeCf1*lunarModelTimeDays;
            
            # Compute the mean elongation of the moon from the sun at the current time.
            lunarElongation = this.ElongCf0 + this.ElongCf1*lunarModelTimeDays;
            
            # Compute the mean anomaly of the sun at the current time.
            solarMeanAnomaly = this.SolAnomCf0 + this.SolAnomCf1*lunarModelTimeDays;
            
            # Compute the inclination of the ecliptic with respect to the mean Earth equator of date, at the current time.  [Actually, this is an error – to calculate the lunar position with respect to the J2000 mean of date coordinate system, we should use the mean obliquity of the ecliptic at the J2000 epoch, meaning the variable ObliquityEcliptic should just be a constant.  However, we can just set OblEclCf0 to the desired constant and OblEclCf1 to zero and achieve the same thing without changing the code.  Since the value of the ecliptic varies in the third or fourth decimal place, this improvement will not make a significant difference.]
            obliquityEcliptic = this.OblEclCf0 + this.OblEclCf1*lunarModelTimeDays;
            
            # Calculate AM
            AM = lunarMeanLongitude - lunarMeanPerigee;
            
            # Brown identified approximately 1600 periodic terms in the motion of the Moon, but only about 10# of these have amplitudes in excess of 1 arcminute.  This routine uses the ten largest solar pertubations and the two largest eccentric terms to correct the mean Longitude of the Moon.
            S1 = lunarMeanLongitude - lunarMeanAscNode;
            perturbation = 6.2887*sin(AM) ...
                         - 1.274*sin(AM - 2*lunarElongation) ...
                         + 0.6583*sin(2*lunarElongation) ...
                         - 0.1855*sin(solarMeanAnomaly) ...
                         - 0.1143*sin(2*S1) ...
                         + 0.053*sin(AM+2*lunarElongation) ...
                         - 0.046*sin(solarMeanAnomaly-2*lunarElongation) ...
                         - 0.035*sin(lunarElongation) ...
                         + 0.213*sin(2*AM) ...
                         - 0.059*sin(2*(AM-lunarElongation) ) ...
                         - 0.057*sin(AM+solarMeanAnomaly - 2*lunarElongation) ...
                         + 0.041*sin(AM-solarMeanAnomaly);
            corLunarMeanLongitude = lunarMeanLongitude + perturbation*this.DEG_TO_RAD;
            
            # Calculate S2.
            S2 = corLunarMeanLongitude - lunarMeanAscNode;
            
            # Calculate the celestial latitude and longitude, lambda and beta, respectively.
            lambda = lunarMeanAscNode + atan2 ( sin(S2) * cos(this.Inclination), cos(S2) );
            beta = atan ( sin(lambda - lunarMeanAscNode)*tan(this.Inclination) );
            
            # Calculate the sine and cosine of the current obliquity of the ecliptic.
            cosObliquity = cos(obliquityEcliptic);
            sinObliquity = sin(obliquityEcliptic);
            
            # Calculate number of centuries since J2000. Calculate precession corrections to Beta and Lambda, the ecliptic longitude and latitude, reducing the coordinates to reference to the J2000 mean of date coordinate system.
            timeJ2000 = (JD - this.JulTime2000) / this.DAYS_PER_CENTURY;
            timeJ2000Sq = timeJ2000*timeJ2000;
            A = this.J2000A1*timeJ2000 + this.J2000A2*timeJ2000Sq;
            B = this.J2000B1*timeJ2000 + this.J2000B2*timeJ2000Sq;
            Cprime = this.J2000Cprim1 + this.J2000Cprim2*timeJ2000 + this.J2000Cprim3*timeJ2000Sq;
            
            betaCor = beta - B*sin(lambda + Cprime);
            lambdaCor = lambda - A + B*cos(lambda+Cprime)*tan(betaCor);
            
            # Calculate parallax (in arc-seconds) based on Brown’s cosine series, including all terms of significance greater than 1 arcsecond.  Determine the Earth to Moon distance.
            parallax = 3422.7 + 186.5398*cos(AM) + 34.3117*cos(AM-2*lunarElongation) ...
                     + 28.2373*cos(2*lunarElongation) + 10.1657*cos(2*AM) ...
                     + 3.0861*cos(AM+2*lunarElongation) ...
                     + 1.9178*cos(solarMeanAnomaly - 2*lunarElongation) ...
                     + 1.4437*cos(AM+solarMeanAnomaly - 2*lunarElongation) ...
                     + 1.1528*cos(AM-solarMeanAnomaly);
            distFromEarthToMoon = this.ErthRad / sin(parallax*this.ARCSEC_TO_RAD); # [km]
            
            # Calculate the unit vector from the center of the Earth to the current moon position in GCI., by calculating the position unit vector in Geocentric Ecliptic Coordinates, then applying a rotation about the Xinertial axis (vernal equinox) by the obliquity of the ecliptic.
            lunarUnitVec    = [0 0 0]';
            lunarUnitVec(1) = cos(betaCor)*cos(lambdaCor);
            lunarUnitVec(2) = cosObliquity*cos(betaCor)*sin(lambdaCor) - sinObliquity*sin(betaCor);
            lunarUnitVec(3) = sinObliquity*cos(betaCor)*sin(lambdaCor) + cosObliquity*sin(betaCor);
            moonPosGci = distFromEarthToMoon * lunarUnitVec; # [km]
            
            # Compute the lunar velocity. Requires two cycles of position and time.
            isValid = false;
            if ~this.isInit(2)
                lunarDeltaTime = moonModelTimeSecs - this.prevMoonModelTimeSecs; # [sec]
                moonVelGci = (moonPosGci - this.prevMoonPosGci) / lunarDeltaTime; # [km/s]
                if ( abs(lunarDeltaTime) > 0)
                    isValid = true;
                end
            else
                moonVelGci = [0 0 0]';
                this.isInit(2) = false;
            end#if            

            # update states
            this.prevMoonPosGci = moonPosGci;
            this.prevMoonModelTimeSecs = moonModelTimeSecs;
        end# lunarModel()
        
        
        function R_TodToJ2000 = TodToJ2000Calculation(this,JD)
            ## True-Of-Date to inertial J2000 frame DCM computation
            T = (JD-this.JulTime2000)/this.DAYS_PER_CENTURY; # Julian century
            
            # Precession matrix (P) calculation
            precessionMtx = eye(3,3);
            zeta  = (2306.2181*T + 0.30188*T^2 + 0.017998*T^3) * this.ARCSEC_TO_RAD; # [rad]
            Z     = (2306.2181*T + 1.09468*T^2 + 0.018203*T^3) * this.ARCSEC_TO_RAD; # [rad]
            theta = (2004.3109*T - 0.42665*T^2 - 0.041833*T^3) * this.ARCSEC_TO_RAD; # [rad]
            precessionMtx(1,1) = -sin(zeta)*sin(Z) + cos(zeta)*cos(Z)*cos(theta);
            precessionMtx(2,1) = -cos(zeta)*sin(Z) - sin(zeta)*cos(Z)*cos(theta);
            precessionMtx(3,1) = -cos(Z)*sin(theta);
            precessionMtx(1,2) = sin(zeta)*cos(Z) + cos(zeta)*sin(Z)*cos(theta);
            precessionMtx(2,2) = cos(zeta)*cos(Z) - sin(zeta)*sin(Z)*cos(theta);
            precessionMtx(3,2) = -sin(Z)*sin(theta);
            precessionMtx(1,3) = cos(zeta)*sin(theta);
            precessionMtx(2,3) = -sin(zeta)*sin(theta);
            precessionMtx(3,3) = cos(theta);
            
            # Nutation matrix (N) calculation
            nutationMtx = eye(3,3);
            Omega  = (125.04452- 1934.136261*T + 0.0020708*T^2 + T^3/450000) * this.DEG_TO_RAD; # [rad]
            L      = (280.4665 + 36000.7698 *T) * this.DEG_TO_RAD; # [rad]
            Lprime = (218.3165 + 481267.8813*T) * this.DEG_TO_RAD; # [rad]
            eps_m  = (23.43929111 - 0.0130047*T - 0.0000001639*T^2 + 0.0000005036*T^3) * this.DEG_TO_RAD; # [rad]
            deps   = (9.20*cos(Omega) + 0.57*cos(2*L) + 0.10*cos(2*Lprime) - 0.09*cos(2*Omega)) * this.ARCSEC_TO_RAD; # [rad]		
			eps_t  = eps_m + deps; #[rad]
            dPsi   = (-17.20*sin(Omega) - 1.32*sin(2*L) - 0.23*sin(2*Lprime) + 0.21*sin(2*Omega)) * this.ARCSEC_TO_RAD; # [rad]           
            nutationMtx(1,1) = cos(dPsi);
            nutationMtx(2,1) = -sin(dPsi)*cos(eps_m);
            nutationMtx(3,1) = -sin(dPsi)*sin(eps_m);
            nutationMtx(1,2) = sin(dPsi)*cos(eps_t);
            nutationMtx(2,2) = cos(eps_t)*cos(eps_m)*cos(dPsi) + sin(eps_t)*sin(eps_m);
            nutationMtx(3,2) = cos(eps_t)*sin(eps_m)*cos(dPsi) - sin(eps_t)*cos(eps_m);
            nutationMtx(1,3) = sin(dPsi)*sin(eps_t);
            nutationMtx(2,3) = sin(eps_t)*cos(eps_m)*cos(dPsi) - cos(eps_t)*sin(eps_m);
            nutationMtx(3,3) = sin(eps_t)*sin(eps_m)*cos(dPsi) + cos(eps_t)*cos(eps_m);          
            
            R_TodToJ2000 = precessionMtx * nutationMtx; # P*N = DCM_TOD2J2000
        end# TodToJ2000Calculation()
        
        
        function [R_EcefToEci, GHA] = EcefToEciCalculation(this, JD, R_TodToJ2000)
            ## ECEF To ECI (J2000) DCM computation 
            # compute Greenwich Mean Sidereal Time (GMST)
            # Ref: https://aa.usno.navy.mil/faq/docs/GAST.php
            # This simple algorithm for computing apparent sidereal time  
            # to an accuracy of ~0.1 sec, equiv to ~1.5 arcsec on the sky. 
            D  = JD - this.JulTime2000;  # JD-2451545.0;
            T  = D/this.DAYS_PER_CENTURY;  # number of centuries since the year 2000
            JD_frac = mod(JD,1);
            if JD_frac>=0.5 # 0h~12h
                JD0 = (JD-JD_frac) + 0.5;
            else # 12h~24h
                JD0 = (JD-JD_frac) - 0.5;
            end
            D0 = JD0 - this.JulTime2000;
            H  = (JD-JD0)*24;   
            #GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T^2; # [hr]
            GMST = this.GmstCf0 + this.GmstCfD0*D0 + this.GmstCfH*H + this.GmstCfT2*T^2; # [hr]
            GHA = mod(GMST,24)*15*this.DEG_TO_RAD; # [rad] 1hr=15deg, Greenwich Hour Angle
            
            # compute earth rotation angle transformation to ECEF
            R_TodToEcef = eye(3);
            R_TodToEcef(1,1) = cos(GHA);
            R_TodToEcef(1,2) = sin(GHA);
            R_TodToEcef(1,3) = 0.0;
            R_TodToEcef(2,1) = -sin(GHA);
            R_TodToEcef(2,2) = cos(GHA);
            R_TodToEcef(2,3) = 0.0;
            R_TodToEcef(3,1) = 0.0;
            R_TodToEcef(3,2) = 0.0;
            R_TodToEcef(3,3) = 1.0;
            
            # Convert from ECEF to ECI(J2000)
            R_EcefToEci = R_TodToJ2000 * R_TodToEcef';
        end# EciToEcefCalculation()
    end#method - protected
        
    ### ----- I/O setting -----
    methods(Access = protected)
        function icon = getIconImpl(~)
            icon = sprintf('');
        end
    end#method - protected

    ### ----- Mask setting -----
    methods(Static, Access = protected)
         function group = getPropertyGroupsImpl(~)
            sunTab = matlab.system.display.SectionGroup(...
                'Title','Solar Model',...
                'PropertyList',{'SunEpoch','CosSunIncl','SinSunIncl','SunMeanCoef1','SunMeanCoef2','SolLongCoef1','SolLongCoef2','SunAnomCrCf1','SunAnomCrCf2','ErthSunDst0','ErthSunDst1','ErthSunDst2','ErthSunDst3'});
            moonTab = matlab.system.display.SectionGroup(...
                'Title','Lunar Model',...
                'PropertyList',{'Inclination','AscNodeCf0','AscNodeCf1','LongCoef0','LongCoef1','PerigeeCf0','PerigeeCf1','ElongCf0','ElongCf1','SolAnomCf0','SolAnomCf1','OblEclCf0','OblEclCf1','ErthRad','JulTimeCoef','JulTime2000','J2000A1','J2000A2','J2000B1','J2000B2','J2000Cprim1','J2000Cprim2','J2000Cprim3'});
            constgp = matlab.system.display.Section(...
                'Title','Constant',...
                'PropertyList',{'SECS_PER_DAY','DAYS_PER_CENTURY','DEG_TO_RAD','ARCSEC_TO_RAD'});
            mainTab = matlab.system.display.SectionGroup(...
                'Title','GMST','Section',[constgp],'PropertyList',{'GmstCf0','GmstCfD0','GmstCfH','GmstCfT2'});             ##ok<NBRAK>
            group = [mainTab,sunTab,moonTab];
         end# getPropertyGroupsImpl()
    end# method - static & protected
end#classdef

### helper functions
function y = wrapInTwoPi(u)
## wrap input angle between 0 <= y < 2*pi
y = u;
while y<0
    y = y + 2*pi;   
end
while y>=2*pi
    y = y - 2*pi;
end
end# wrapInTwoPi()
