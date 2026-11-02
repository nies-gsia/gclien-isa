/*
 * The MIT License (MIT)
 * 
 * Copyright (c) 2016-2017 Max Lieblich
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *   
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

/*
 * For the original file, see "wbgt.c.original" in the src directory of this 
 * package. That file has been modified in the following ways:
 * 
 * 1. The "main" function has been removed.
 * 2. All function declarations have been converted from K&R style to ANSI C99
 *    compliant style and moved to "wbgt.h".
 * 3. All floats have been converted to doubles to integrate more easily with R.
 * 4. Parameter documentation for functions is included in "wbgt.h"; these are
 *    copies of comments in "wbgt.c.original" accompanying the original K&R
 *    parameter declarations.
 * 5. For documentation in its original location, see "wbgt.c.original"
 * 
 * Max Lieblich, University of Washington, 2016
 */

/* Original license follows */

/*
               Copyright © 2008, UChicago Argonne, LLC
                       All Rights Reserved

                        WBGT, Version 1.1

			     James C. Liljegren
              Decision & Information Sciences Division

			     OPEN SOURCE LICENSE

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.  Software changes, 
   modifications, or derivative works, should be noted with comments and 
   the author and organization’s name.

2. Redistributions in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the documentation 
   and/or other materials provided with the distribution.

3. Neither the names of UChicago Argonne, LLC or the Department of Energy 
   nor the names of its contributors may be used to endorse or promote products 
   derived from this software without specific prior written permission.

4. The software and the end-user documentation included with the 
   redistribution, if any, must include the following acknowledgment:

   "This product includes software produced by UChicago Argonne, LLC 
   under Contract No. DE-AC02-06CH11357 with the Department of Energy.”

******************************************************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, 
NOR UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS 
OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, 
COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS, PRODUCT, OR 
PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
*******************************************************************************************
*/

/*
 * Modifications copyright (c) 2026 Jun'ya Takakura,
 * National Institute for Environmental Studies (NIES).
 *
 * The modifications made by NIES are licensed under the BSD 3-Clause License
 * (SPDX-License-Identifier: BSD-3-Clause; see the LICENSE file). The original
 * code remains under the licenses above.
 *
 * Summary of the NIES modifications (individual changes are noted in the
 * comments below as "Updated by NIES" or "Author: NIES"):
 *  - The wet-bulb temperature is solved with the Newton-Raphson method
 *    (new helper functions such as Twb_rhs).
 *  - The psychrometric wet-bulb temperature is not calculated when the
 *    argument Tpsy is NULL.
 *  - Ported from C to C++ (modernized C-style code; #define constants and
 *    max/min macros replaced with constexpr and std::max/std::min).
 */
 
/*
 *  Purpose: to demonstrate the use of the subroutine calc_wbgt to calculate
 *           the wet bulb-globe temperature (WBGT).  The program reads input 
 *           data from a file containing meteorological measurements then 
 *           calls calc_wbgt to compute the WBGT.
 *
 *           The inputs and outputs are fully described in calc_wbgt.
 *
 *  Author:  James C. Liljegren
 *		 Decision and Information Sciences Division
 *		 Argonne National Laboratory
 *
 *  Updated by NIES:
 *       When th argument Tspy is NULL, Tspy is not calculated.
 */		
 
#include	<cmath>
#include	<algorithm>
#include	<array>
/* 
 * header added to factor out the calc_wbgt function
 * Max Lieblich, University of Washington, 2016
 */
#include "wbgtn.h"

static double Twb_rhs(double Tw, double Tair, double Pair, double speed,
					  int rad, double eair, double Qr, double Fatm2);
static double d_Twb_rhs(double Tw, double Tair, double Pair, double speed, int rad,
						double eair, double Qr, double Fatm2);
static double Tglobe_rhs(double Tglobe_prev, double Tair, double Pair,
						 double speed, double Qr, double term3);
static double d_Tglobe_rhs(double Tglobe_prev, double Tair, double Pair,
						   double speed, double Qr, double term3);

/* ============================================================================
 *  Purpose: to calculate the outdoor wet bulb-globe temperature, which is 
 *           the weighted sum of the air temperature (dry bulb), the globe temperature, 
 *           and the natural wet bulb temperature: Twbg = 0.1 * Tair + 0.7 * Tnwb + 0.2 * Tg.
 *           The program predicts Tnwb and Tg using meteorological input data then combines
 *           the results to produce Twbg.
 *
 *		 Modified 2-Nov-2009: calc_wbgt returns -1 if either subroutines Tg or Tnwb return -9999,
 *		 which signals a failure to converge, probably due to a bad input value; otherwise, calc_wbgt
 *		 returns 0.
 *
 *           If the 2-m wind speed is not available, it is estimated using a wind speed at another level.
 *	
 *  Reference: Liljegren, J. C., R. A. Carhart, P. Lawday, S. Tschopp, and R. Sharp:
 *             Modeling the Wet Bulb Globe Temperature Using Standard Meteorological
 *             Measurements. The Journal of Occupational and Environmental Hygiene, 
 *             vol. 5:10, pp. 645-655, 2008.
 *
 *  Author:  James C. Liljegren 
 *		 Decision and Information Sciences Division
 *		 Argonne National Laboratory
 */


int calc_wbgt(int year, int month, int day, int hour, int minute, int gmt, int avg, 
              double lat, double lon, double solar, double pres, double Tair, double relhum, 
              double speed, double zspeed, double dT, int urban, double *est_speed,
			  double *Tg, double *Tnwb, double *Tpsy, double *Twbg)
{
	double	cza,	/* cosine of solar zenith angle						*/
		fdir,	/* fraction of solar irradiance due to direct beam			*/
		tk,	/* temperature converted to kelvin						*/
		rh;	/* relative humidity, fraction between 0 and 1				*/

	double hour_gmt, dday;
	
	bool daytime;
	int stability_class; 

	/* 
 *  convert time to GMT and center in avg period;
 */
	hour_gmt = hour - gmt + ( minute - 0.5 * avg ) / 60.; 
	dday = day + hour_gmt / 24.;
/*
 *  calculate the cosine of the solar zenith angle and fraction of solar irradiance
 *  due to the direct beam; adjust the solar irradiance if it is out of bounds
 */
	calc_solar_parameters(year, month, dday, lat, lon, &solar, &cza, &fdir);	
/* 
 *  estimate the wind speed, if necessary
 */
	if ( zspeed != REF_HEIGHT ) {	
		daytime = (cza > 0.);
		stability_class = stab_srdt(daytime, speed, solar, dT);
		*est_speed = est_wind_speed(speed, zspeed, stability_class, urban);
		speed = *est_speed;
	}
/*
 *  unit conversions
 */
	tk = Tair + 273.15; /* degC to kelvin */
	rh = 0.01 * relhum; /* % to fraction  */
/*
 *  calculate the globe, natural wet bulb, psychrometric wet bulb, and 
 *  outdoor wet bulb globe temperatures
 */
	*Tg   = Tglobe(tk, rh, pres, speed, solar, fdir, cza);
	if ( *Tg == -9999 ) {
		*Twbg = -9999;
		return -1;
	}
	*Tnwb = Twb(tk, rh, pres, speed, solar, fdir, cza, 1);
	if ( *Tnwb == -9999 ) {
		*Twbg = -9999;
		return -1;
	}
	if ( Tpsy ) *Tpsy = Twb(tk, rh, pres, speed, solar, fdir, cza, 0);
	*Twbg = 0.1 * Tair + 0.2 * (*Tg) + 0.7 * (*Tnwb); 
		
	return 0;
}

/* ============================================================================
 *  Purpose: to calculate the cosine solar zenith angle and the fraction of the
 *		 solar irradiance due to the direct beam.
 *	
 *  Author:  James C. Liljegren 
 *		 Decision and Information Sciences Division
 *		 Argonne National Laboratory
 */

int	calc_solar_parameters(int year, int month, double day, double lat, double lon, 
                          double *solar, double *cza, double *fdir)
{
	double	toasolar, normsolar; 
	
	double days_1900 = 0.0, ap_ra, ap_dec, elev, refr, azim, soldist;
	
	solarposition(year, month, day, days_1900, lat, lon,
		&ap_ra, &ap_dec, &elev, &refr, &azim, &soldist);
	*cza = cos( (90.-elev)*DEG_RAD );
	toasolar = SOLAR_CONST * std::max(0.,*cza) / (soldist*soldist);
/* 
 *  if the sun is not fully above the horizon 
 *  set the maximum (top of atmosphere) solar = 0
 */	
	if ( *cza < CZA_MIN ) toasolar = 0.;
	if ( toasolar > 0. ) {
/*
 *  account for any solar sensor calibration errors and
 *  make the solar irradiance consistent with normsolar
 */
		normsolar = std::min( *solar/toasolar, NORMSOLAR_MAX );
		*solar = normsolar * toasolar;
/*
 *  calculate the fraction of the solar irradiance due to the direct beam
 */
		if ( normsolar > 0. ) {
			*fdir = exp( 3. - 1.34 * normsolar - 1.65 / normsolar );
			*fdir = std::max(std::min(*fdir, 0.9), 0.0);
		} 
		else
			*fdir = 0.;
	}
	else
		*fdir = 0.;

	return 0;
}


int	calc_solar_cza(int year, int month, double day, double lat, double lon, double* cza)
{
	double days_1900 = 0.0, ap_ra, ap_dec, elev, refr, azim, soldist;

	solarposition(year, month, day, days_1900, (double)lat, (double)lon,
		&ap_ra, &ap_dec, &elev, &refr, &azim, &soldist);
	*cza = cos((90. - elev) * DEG_RAD);

	return 0;
}

/* ============================================================================
 *  Purpose: to calculate the natural wet bulb temperature.
 *	
 *  Author:  James C. Liljegren 
 *		 Decision and Information Sciences Division
 *		 Argonne National Laboratory
 *
 *  Updated by NIES:
 *       Newton-Raphson method is used as an iterative solution.
 */
 
double Twb(double Tair, double rh, double Pair, double speed, 
		   double solar, double fdir, double cza, int rad)
{
	double sza, Tsfc, Tdew, Tw, eair, tan_sza;
	double Qr; /* thermal radiation emitted by the atmosphere and ground surface
				  in units of STEFANB */
	double Fatm2; /* 2nd term of Eq.(12) */
	double dT, f, df;
	bool converged;
	int iter;
	
	Tsfc = Tair;
	sza = acos(cza); /* solar zenith angle, radians */
	eair = rh * esat(Tair,0);
	Tdew = dew_point(eair,0);
	Tw = Tdew; /* first guess is the dew point temperature */
	Qr = EMIS_SFC * pow(Tsfc,4.) + emis_atm(Tair, rh) * pow(Tair,4.);
	tan_sza = tan(sza);
	Fatm2 = (1.-ALB_WICK) * solar *
		( (1.-fdir)*(1.+0.25*D_WICK/L_WICK) + fdir*((tan_sza/PI)+0.25*D_WICK/L_WICK) + ALB_SFC );
	converged = false;
	for (iter = 1; iter <= MAX_ITER; iter++) {
		f = Tw - Twb_rhs(Tw, Tair, Pair, speed, rad, eair, Qr, Fatm2);
		df = 1 - d_Twb_rhs(Tw, Tair, Pair, speed, rad, eair, Qr, Fatm2);
		dT = - f / df;
		Tw += dT;
		if ( std::abs(dT) < CONVERGENCE ) {
			converged = true;
			break;
		}
	}	

	if ( converged ) 
		return Tw - 273.15;
	else
		return -9999.;
}

/* ============================================================================
 * Purpose: to calculate the convective heat transfer coefficient in W/(m2 K)
 *          for a long cylinder in cross flow.
 *
 * Reference: Bedingfield and Drew, eqn 32 
 *
 * Updated by NIES
 *      The constant value (PowPr) is added as an static const variable.
 */
 
double h_cylinder_in_air(double diameter, double length, double Tair, double Pair, double speed)
{
	static const double a = 0.56,  /* parameters from Bedingfield and Drew */
			 b = 0.281,
			 c = 0.4;
	static const double PowPr = pow(Pr,(1.-a));
	
	double	density,
		Re,	/* Reynolds number								*/
		Nu;	/* Nusselt number									*/

	density = Pair * 100. / ( R_AIR * Tair );
	Re = std::max(speed,MIN_SPEED) * density * diameter / viscosity(Tair);
	Nu = b * pow(Re,(1.-c)) * PowPr;
	return Nu * thermal_cond(Tair) / diameter;
}

/* ============================================================================
 *  Purpose: to calculate the globe temperature.
 *	
 *  Author:  James C. Liljegren 
 *		 Decision and Information Sciences Division
 *		 Argonne National Laboratory
 *
 *  Updated by NIES:
 *       Newton-Raphson method is used as an iterative solution.
 */
 
double Tglobe(double Tair, double rh, double Pair, double speed, double solar, double fdir, double cza)
{
	double Tsfc, Tg;
	double Qr; /* thermal radiation emitted by the atmosphere and ground surface
				  in units of STEFANB */
	double term3; /* the third term of the right-hand side of Eq.(17) */
	double dT, f, df, Tg3;
	bool converged;
	int iter;
	
	Tsfc = Tair;
	Tg = Tair; /* first guess is the air temperature */
	Qr = emis_atm(Tair,rh)*pow(Tair,4.) + EMIS_SFC*pow(Tsfc,4.);
	term3 = solar/(2.*STEFANB*EMIS_GLOBE)*(1.-ALB_GLOBE)*(fdir*(1./(2.*cza)-1.)+1.+ ALB_SFC);	
	converged = false;
	for (iter = 1; iter <= MAX_ITER; iter++) {
		Tg3 = Tg * Tg * Tg;
		f = Tg3 * Tg;
		f -= Tglobe_rhs(Tg, Tair, Pair, speed, Qr, term3);
		df = 4 * Tg3 - d_Tglobe_rhs(Tg, Tair, Pair, speed, Qr, term3);
		dT = - f / df;
		Tg += dT;
		if ( std::abs(dT) < CONVERGENCE ) {
			converged = true;
			break;
		}
	}
	
	if ( converged ) 
		return Tg - 273.15;
	else
		return -9999.;
}

/* ============================================================================
 * Purpose: to calculate the convective heat transfer coefficient, W/(m2 K)
 *          for flow around a sphere.
 *
 * Reference: Bird, Stewart, and Lightfoot (BSL), page 409. 
 *
 * Updated by NIES
 *      The constant value (PowPr) is added as an static const variable.
 */
 
double h_sphere_in_air(double diameter, double Tair, double Pair, double speed)
{
	static const double PowPr = pow(Pr,0.3333);
	double	density,
		Re,	/* Reynolds number							*/
		Nu;	/* Nusselt number								*/
	
	density = Pair * 100. / ( R_AIR * Tair );
	Re = std::max(speed,MIN_SPEED) * density * diameter / viscosity(Tair);
	Nu = 2.0 + 0.6 * sqrt(Re) * PowPr;
	return Nu * thermal_cond(Tair) / diameter;
}
 

/* ============================================================================
 *  Purpose: calculate the saturation vapor pressure (mb) over liquid water
 *           (phase = 0) or ice (phase = 1).
 *
 *  Reference: Buck's (1981) approximation (eqn 3) of Wexler's (1976) formulae.
 */
 
double esat(double tk, int phase)
{
	double y, es;
	
	if ( phase == 0 ) {	/* over liquid water */
		y = (tk - 273.15)/(tk - 32.18);
		es = 6.1121 * exp( 17.502 * y );
//		es = (1.0007 + (3.46E-6 * pres)) * es /* correction for moist air, if pressure is available */
	} 
	else {			/* over ice */
		y = (tk - 273.15)/(tk - 0.6);
		es = 6.1115 * exp( 22.452 * y );
//		es = (1.0003 + (4.18E-6 * pres)) * es /* correction for moist air, if pressure is available */
	}
	
	es = 1.004 * es;  /* correction for moist air, if pressure is not available; for pressure > 800 mb */
//	es = 1.0034 * es; /* correction for moist air, if pressure is not available; for pressure down to 200 mb */

	return es;
}

/* ============================================================================
 *  Purpose: calculate the dew point (phase=0) or frost point (phase=1) 
 *           temperature, K.
 */
 
double dew_point(double e, int phase)
{
	double z, tdk;
	
	if ( phase == 0 ) {	/* dew point */
		z = log( e / (6.1121*1.004) );
		tdk = 273.15 + 240.97*z/(17.502-z);
	}
	else {	            /* frost point */
		z = log( e / (6.1115*1.004) );
		tdk = 273.15 + 272.55*z/(22.452-z);
	}
	
	return tdk;
}
 
/* ============================================================================
 *  Purpose: calculate the viscosity of air, kg/(m s)
 *
 *  Reference: BSL, page 23.
 */
 
double viscosity(double Tair)
{
	static const double sigma = 3.617,
		     eps_kappa = 97.0;
			 
	double	Tr, omega;
	
	Tr = Tair / eps_kappa;
	omega = ( Tr - 2.9 ) / 0.4 * ( -0.034 ) + 1.048;
	return 2.6693E-6 * sqrt( M_AIR*Tair ) / ( sigma * sigma * omega );
}

/* ============================================================================
 *  Purpose: calculate the thermal conductivity of air, W/(m K)
 *
 *  Reference: BSL, page 257.
 */
 
double thermal_cond(double Tair)
{			 
	return ( Cp + 1.25 * R_AIR ) * viscosity(Tair);
}

/* ============================================================================
 *  Purpose: calculate the diffusivity of water vapor in air, m2/s
 *
 *  Reference: BSL, page 505.
 *
 *  Updated by NIES
 *      The values (Pcrit13, Tcrit512, Tcrit12, Mmix) are changed from local
 *      variables to static const ones.
 */
 
double diffusivity(double Tair, double Pair)
{
	static const double Pcrit_air = 36.4, 
			 Pcrit_h2o = 218., 
			 Tcrit_air = 132., 
			 Tcrit_h2o = 647.3,
			 a = 3.640E-4, 
			 b = 2.334;			 
	static const double Pcrit13 = pow( ( Pcrit_air * Pcrit_h2o ),(1./3.) );
	static const double Tcrit512 = pow( ( Tcrit_air * Tcrit_h2o ),(5./12.) );
	static const double Tcrit12 = sqrt( Tcrit_air * Tcrit_h2o );
	static const double Mmix = sqrt( 1./M_AIR + 1./M_H2O );
	double Patm;

	Patm = Pair / 1013.25 ; /* convert pressure from mb to atmospheres */
	
	return a * pow( (Tair/Tcrit12),b) * Pcrit13 * Tcrit512 * Mmix / Patm * 1E-4;
}

/* ============================================================================
 *  Purpose: calculate the heat of evaporation, J/(kg K), for temperature
 *           in the range 283-313 K.
 *
 *  Reference: Van Wylen and Sonntag, Table A.1.1
 */
 
double evap(double Tair)
{			 
	return (313.15 - Tair)/30. * (-71100.) + 2.4073E6;
}

/* ============================================================================
 *  Purpose: calculate the atmospheric emissivity.
 *
 *  Reference: Oke (2nd edition), page 373.
 */
 
double emis_atm(double Tair, double rh)
{
	double e; 

	e = rh * esat(Tair,0);
	return 0.575 * pow(e, 0.143);
}

/* ============================================================================
 *  Version 3.0 - February 20, 1992.
 *
 *  solarposition() employs the low precision formulas for the Sun's coordinates
 *  given in the "Astronomical Almanac" of 1990 to compute the Sun's apparent
 *  right ascension, apparent declination, altitude, atmospheric refraction
 *  correction applicable to the altitude, azimuth, and distance from Earth.
 *  The "Astronomical Almanac" (A. A.) states a precision of 0.01 degree for the
 *  apparent coordinates between the years 1950 and 2050, and an accuracy of
 *  0.1 arc minute for refraction at altitudes of at least 15 degrees.
 *
 *  The following assumptions and simplifications are made:
 *  -> refraction is calculated for standard atmosphere pressure and temperature
 *     at sea level.
 *  -> diurnal parallax is ignored, resulting in 0 to 9 arc seconds error in
 *     apparent position.
 *  -> diurnal aberration is also ignored, resulting in 0 to 0.02 second error
 *     in right ascension and 0 to 0.3 arc second error in declination.
 *  -> geodetic site coordinates are used, without correction for polar motion
 *     (maximum amplitude of 0.3 arc second) and local gravity anomalies.
 *  -> local mean sidereal time is substituted for local apparent sidereal time
 *     in computing the local hour angle of the Sun, resulting in an error of
 *     about 0 to 1 second of time as determined explicitly by the equation of
 *     the equinoxes.
 *
 *  Right ascension is measured in hours from 0 to 24, and declination in
 *  degrees from 90 to -90.
 *  Altitude is measured from 0 degrees at the horizon to 90 at the zenith or
 *  -90 at the nadir. Azimuth is measured from 0 to 360 degrees starting at
 *  north and increasing toward the east at 90.
 *  The refraction correction should be added to the altitude if Earth's
 *  atmosphere is to be accounted for.
 *  Solar distance from Earth is in astronomical units, 1 a.u. representing the
 *  mean value.
 *
 *  The necessary input parameters are:
 *  -> the date, specified in one of three ways:
 *       1) year, month, day.fraction
 *       2) year, daynumber.fraction
 *       3) days.fraction elapsed since January 0, 1900.
 *  -> site geodetic (geographic) latitude and longitude.
 *
 *  Refer to the function declaration for the parameter type specifications and
 *  formats.
 *
 *  solarposition() returns -1 if an input parameter is out of bounds, or 0 if
 *  values were written to the locations specified by the output parameters.
 */

/*  Author: Nels Larson
 *          Pacific Northwest National Laboratory
 *          P.O. Box 999
 *          Richland, WA 99352
 *          U.S.A.
 */


int solarposition(int year, int month, double day, double days_1900, double latitude, 
                  double longitude, double *ap_ra, double *ap_dec, double *altitude, 
                  double *refraction, double *azimuth, double *distance)
{
  /*
   * Remove unnecessary declaration; moved to wgbt.h
   * Max Lieblich, University of Washington, 2016
   */
  /*int    daynum();        Computes a sequential daynumber during a year. */

  int    daynumber,       /* Sequential daynumber during a year. */
         delta_days,      /* Whole days since 2000 January 0. */
         delta_years;     /* Whole years since 2000. */
  double cent_J2000,      /* Julian centuries since epoch J2000.0 at 0h UT. */
         cos_alt,         /* Cosine of the altitude of Sun. */
         cos_apdec,       /* Cosine of the apparent declination of Sun. */
         cos_az,          /* Cosine of the azimuth of Sun. */
         cos_lat,         /* Cosine of the site latitude. */
         cos_lha,         /* Cosine of the local apparent hour angle of Sun. */
         days_J2000,      /* Days since epoch J2000.0. */
         ecliptic_long,   /* Solar ecliptic longitude. */
         lmst,            /* Local mean sidereal time. */
         local_ha,        /* Local mean hour angle of Sun. */
         gmst0h,          /* Greenwich mean sidereal time at 0 hours UT. */
         integral,        /* Integral portion of double precision number. */
         mean_anomaly,    /* Earth mean anomaly. */
         mean_longitude,  /* Solar mean longitude. */
         mean_obliquity,  /* Mean obliquity of the ecliptic. */
         pressure =       /* Earth mean atmospheric pressure at sea level */
           1013.25,       /*   in millibars. */
         sin_apdec,       /* Sine of the apparent declination of Sun. */
         sin_az,          /* Sine of the azimuth of Sun. */
         sin_lat,         /* Sine of the site latitude. */
         tan_alt,         /* Tangent of the altitude of Sun. */
         temp =           /* Earth mean atmospheric temperature at sea level */
           15.0,          /*   in degrees Celsius. */
         ut;              /* UT hours since midnight. */


  /* Check latitude and longitude for proper range before calculating dates.
   */
  if (latitude < -90.0 || latitude > 90.0 ||
      longitude < -180.0 || longitude > 180.0)
    return -1;

  /* If year is not zero then assume date is specified by year, month, day.
   * If year is zero then assume date is specified by days_1900.
   */
  if (year != 0)
  /* Date given by {year, month, day} or {year, 0, daynumber}. */
  {
    /* The low-precision Astronomical Almanac formulas are stated valid for
     * 1950-2049. For dates outside this range (e.g. ISIMIP3b historical/SSP
     * runs covering 1850-2100) the Sun's apparent position on a given calendar
     * day differs only negligibly from a nearby in-range year -- the dominant
     * variation is the leap-year phase (< ~0.5 deg in declination) -- so clamp
     * the year used for the geometry instead of failing. The caller's real
     * year is unaffected (it is passed by value), so output naming and bias
     * correction still use the true year. */
    if (year < 1950) year = 1950;
    else if (year > 2049) year = 2049;
    if (month != 0)
    {
      if (month < 1 || month > 12 || day < 0.0 || day > 33.0)
        return -1;

      daynumber = daynum(year, month, (int)day);
    }
    else
    {
      if (day < 0.0 || day > 368.0)
        return -1;

      daynumber = (int)day;
    }

    /* Construct Julian centuries since J2000 at 0 hours UT of date,
     * days.fraction since J2000, and UT hours.
     */
    delta_years = year - 2000;
    /* delta_days is days from 2000/01/00 (1900's are negative). */
    delta_days = delta_years * 365 + delta_years / 4 + daynumber;
    if (year > 2000)
      delta_days += 1;
    /* J2000 is 2000/01/01.5 */
    days_J2000 = delta_days - 1.5;

    cent_J2000 = days_J2000 / 36525.0;

    ut = modf(day, &integral);
    days_J2000 += ut;
    ut *= 24.0;
  }
  else
  /* Date given by days_1900. */
  {
    /* days_1900 is 18262 for 1950/01/00, and 54788 for 2049/12/32.
     * A. A. 1990, K2-K4. */
    if (days_1900 < 18262.0 || days_1900 > 54788.0)
      return -1;

    /* Construct days.fraction since J2000, UT hours, and
     * Julian centuries since J2000 at 0 hours UT of date.
     */
    /* days_1900 is 36524 for 2000/01/00. J2000 is 2000/01/01.5 */
    days_J2000 = days_1900 - 36525.5;

    ut = modf(days_1900, &integral) * 24.0;

    cent_J2000 = (integral - 36525.5) / 36525.0;
  }


  /* Compute solar position parameters.
   * A. A. 1990, C24.
   */

  mean_anomaly = (357.528 + 0.9856003 * days_J2000);
  mean_longitude = (280.460 + 0.9856474 * days_J2000);

  /* Put mean_anomaly and mean_longitude in the range 0 -> 2 pi. */
  mean_anomaly = modf(mean_anomaly / 360.0, &integral) * TWOPI;
  mean_longitude = modf(mean_longitude / 360.0, &integral) * TWOPI;

  mean_obliquity = (23.439 - 4.0e-7 * days_J2000) * DEG_RAD;
  ecliptic_long = ((1.915 * sin(mean_anomaly)) +
                   (0.020 * sin(2.0 * mean_anomaly))) * DEG_RAD +
                  mean_longitude;

  *distance = 1.00014 - 0.01671 * cos(mean_anomaly) -
              0.00014 * cos(2.0 * mean_anomaly);

  /* Tangent of ecliptic_long separated into sine and cosine parts for ap_ra. */
  *ap_ra = atan2(cos(mean_obliquity) * sin(ecliptic_long), cos(ecliptic_long));

  /* Change range of ap_ra from -pi -> pi to 0 -> 2 pi. */
  if (*ap_ra < 0.0)
    *ap_ra += TWOPI;
  /* Put ap_ra in the range 0 -> 24 hours. */
  *ap_ra = modf(*ap_ra / TWOPI, &integral) * 24.0;

  *ap_dec = asin(sin(mean_obliquity) * sin(ecliptic_long));

  /* Calculate local mean sidereal time.
   * A. A. 1990, B6-B7.
   */

  /* Horner's method of polynomial exponent expansion used for gmst0h. */
  gmst0h = 24110.54841 + cent_J2000 * (8640184.812866 + cent_J2000 *
            (0.093104 - cent_J2000 * 6.2e-6));
  /* Convert gmst0h from seconds to hours and put in the range 0 -> 24. */
  gmst0h = modf(gmst0h / 3600.0 / 24.0, &integral) * 24.0;
  if (gmst0h < 0.0)
    gmst0h += 24.0;

  /* Ratio of lengths of mean solar day to mean sidereal day is 1.00273790934
   * in 1990. Change in sidereal day length is < 0.001 second over a century.
   * A. A. 1990, B6.
   */
  lmst = gmst0h + (ut * 1.00273790934) + longitude / 15.0;
  /* Put lmst in the range 0 -> 24 hours. */
  lmst = modf(lmst / 24.0, &integral) * 24.0;
  if (lmst < 0.0)
    lmst += 24.0;

  /* Calculate local hour angle, altitude, azimuth, and refraction correction.
   * A. A. 1990, B61-B62.
   */

  local_ha = lmst - *ap_ra;
  /* Put hour angle in the range -12 to 12 hours. */
  if (local_ha < -12.0)
    local_ha += 24.0;
  else if (local_ha > 12.0)
    local_ha -= 24.0;

  /* Convert latitude and local_ha to radians. */
  latitude *= DEG_RAD;
  local_ha = local_ha / 24.0 * TWOPI;

  cos_apdec = cos(*ap_dec);
  sin_apdec = sin(*ap_dec);
  cos_lat = cos(latitude);
  sin_lat = sin(latitude);
  cos_lha = cos(local_ha);

  *altitude = asin(sin_apdec * sin_lat + cos_apdec * cos_lha * cos_lat);

  cos_alt = cos(*altitude);
  /* Avoid tangent overflow at altitudes of +-90 degrees.
   * 1.57079615 radians is equal to 89.99999 degrees.
   */
  if (std::abs(*altitude) < 1.57079615)
    tan_alt = tan(*altitude);
  else
    tan_alt = 6.0e6;

  cos_az = (sin_apdec * cos_lat - cos_apdec * cos_lha * sin_lat) / cos_alt;
  sin_az = -(cos_apdec * sin(local_ha) / cos_alt);
  *azimuth = acos(cos_az);

  /* Change range of azimuth from 0 -> pi to 0 -> 2 pi. */
  if (atan2(sin_az, cos_az) < 0.0)
    *azimuth = TWOPI - *azimuth;

  /* Convert ap_dec, altitude, and azimuth to degrees. */
  *ap_dec *= RAD_DEG;
  *altitude *= RAD_DEG;
  *azimuth *= RAD_DEG;

  /* Compute refraction correction to be added to altitude to obtain actual
   * position.
   * Refraction calculated for altitudes of -1 degree or more allows for a
   * pressure of 1040 mb and temperature of -22 C. Lower pressure and higher
   * temperature combinations yield less than 1 degree refraction.
   * NOTE:
   * The two equations listed in the A. A. have a crossover altitude of
   * 19.225 degrees at standard temperature and pressure. This crossover point
   * is used instead of 15 degrees altitude so that refraction is smooth over
   * the entire range of altitudes. The maximum residual error introduced by
   * this smoothing is 3.6 arc seconds at 15 degrees. Temperature or pressure
   * other than standard will shift the crossover altitude and change the error.
   */
  if (*altitude < -1.0 || tan_alt == 6.0e6)
    *refraction = 0.0;
  else
  {
    if (*altitude < 19.225)
    {
      *refraction = (0.1594 + (*altitude) * (0.0196 + 0.00002 * (*altitude))) *
                    pressure;
      *refraction /= (1.0 + (*altitude) * (0.505 + 0.0845 * (*altitude))) *
                     (273.0 + temp);
    }
    else
      *refraction = 0.00452 * (pressure / (273.0 + temp)) / tan_alt;
  }
/* 
 *  to match Michalsky's sunae program, the following line was inserted
 *  by JC Liljegren to add the refraction correction to the solar altitude
 */
  *altitude = *altitude + *refraction;
  return 0;
}


/* ============================================================================
 * 'daynum()' returns the sequential daynumber of a calendar date during a
 *  Gregorian calendar year (for years 1 onward).
 *  The integer arguments are the four-digit year, the month number, and
 *  the day of month number.
 *  (Jan. 1 = 01/01 = 001; Dec. 31 = 12/31 = 365 or 366.)
 *  A value of -1 is returned if the year is out of bounds.
 */

/* Author: Nels Larson
 *         Pacific Northwest Lab.
 *         P.O. Box 999
 *         Richland, WA 99352
 *         U.S.A.
 */

int daynum(int year, int month, int day)
{
  static constexpr std::array<int, 13> begmonth = {0,0,31,59,90,120,151,181,212,243,273,304,334};
  int dnum;
  bool leapyr = false;

  /* There is no year 0 in the Gregorian calendar and the leap year cycle
   * changes for earlier years. */

  if (year < 1)
    return -1;

  /* Leap years are divisible by 4, except for centurial years not divisible
   * by 400. */

  leapyr = ((year%4) == 0 && (year%100) != 0) || (year%400) == 0;

  dnum = begmonth[month] + day;
  if (leapyr && (month > 2))
    dnum += 1;

  return dnum;
}

/* ============================================================================
 *  Purpose: estimate 2-m wind speed for all stability conditions
 *
 *  Reference: EPA-454/5-99-005, 2000, section 6.2.5
 */

double est_wind_speed(double speed, double zspeed, int stability_class, int urban)
{
	constexpr std::array<double, 6> urban_exp = { 0.15, 0.15, 0.20, 0.25, 0.30, 0.30 };
	constexpr std::array<double, 6> rural_exp = { 0.07, 0.07, 0.10, 0.15, 0.35, 0.55 };
	double exponent, est_speed;
		
	if ( urban )
		exponent = urban_exp[stability_class-1];
	else
		exponent = rural_exp[stability_class-1];
	
	est_speed = speed * pow( REF_HEIGHT/zspeed, exponent );
	est_speed = std::max( est_speed, MIN_SPEED );
	return est_speed;
}
	
/* ============================================================================
 *  Purpose: estimate the stability class
 *
 *  Reference: EPA-454/5-99-005, 2000, section 6.2.5
 */
int stab_srdt(int daytime, double speed, double solar, double dT)
{
	static constexpr std::array<std::array<int, 8>, 6> lsrdt = {{
			{1, 1, 2, 4, 0, 5, 6, 0},
			{1, 2, 3, 4, 0, 5, 6, 0},
			{2, 2, 3, 4, 0, 4, 4, 0},
			{3, 3, 4, 4, 0, 0, 0, 0},
			{3, 4, 4, 4, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0}
	}};
			
	int	i,j;
	
	if ( daytime ) {
		if ( solar >= 925.0 )
			j = 0;
		else if ( solar >= 675.0 )
			j = 1;
		else if ( solar >= 175.0 )
			j = 2;
		else
			j = 3;
			
		if ( speed >= 6.0 ) 
			i = 4;
		else if ( speed >= 5.0 )
			i = 3;
		else if ( speed >= 3.0 )
			i = 2;
		else if ( speed >= 2.0 )
			i = 1;
		else
			i = 0;
	} 
	else {
		if ( dT >= 0.0 )
			j = 6;
		else
			j = 5;
			
		if ( speed >= 2.5 )
			i = 2;
		else if ( speed >= 2.0 )
			i = 1;
		else
			i = 0;
	}
	return lsrdt[i][j];
}

/* ============================================================================
 *  Purpose: to calculate the right-hand side of Eq.(9).
 *	
 *  Author:  NIES
 */

static double Twb_rhs(double Tw, double Tair, double Pair, double speed,
					  int rad, double eair, double Qr, double Fatm2)
{
	static const double a = 0.56; /* from Bedingfield and Drew */	
	double Tref, ewick, density, v;
	double h; 	 /* convective heat transfer coefficient */
	double Fatm; /* radiative heating term */
	double Sc;	 /* Schmidt number */
	
	Tref = 0.5*( Tw + Tair );	/* evaluate properties at the average temperature */
	if (!rad) {
		h = 0.0;
		Fatm = 0.0;
	} else {
		h = h_cylinder_in_air(D_WICK, L_WICK, Tref, Pair, speed);
		Fatm = STEFANB * EMIS_WICK *
			   ( 0.5*Qr - pow(Tw,4.) )
			 + Fatm2;
	}
	ewick = esat(Tw,0);
	density = Pair * 100. / (R_AIR * Tref);
	Sc = viscosity(Tref)/(density*diffusivity(Tref,Pair));

	v = Tair - evap(Tref)/RATIO * (ewick-eair)/(Pair-ewick) * pow(Pr/Sc,a);
	if (rad != 0) v += Fatm/h * rad;

	return v;
}

/* ============================================================================
 *  Purpose: to calculate the first derivative of the right-hand side of Eq.(9).
 *	
 *  Author:  NIES
 */

static double d_Twb_rhs(double Tw, double Tair, double Pair, double speed, int rad,
						double eair, double Qr, double Fatm2)
{
	double v;

	v = Twb_rhs(Tw + DELTA_DERIV, Tair, Pair, speed, rad, eair,
				Qr, Fatm2);
	v -= Twb_rhs(Tw - DELTA_DERIV, Tair, Pair, speed, rad, eair,
				 Qr, Fatm2);
	v /= (2 * DELTA_DERIV);
	
	return v;
}

/* ============================================================================
 *  Purpose: to calculate the right-hand side of Eq.(17).
 *	
 *  Author:  NIES
 */

static double Tglobe_rhs(double Tg, double Tair, double Pair, double speed,
						 double Qr, double term3)
{
	double Tref, h, v;

	Tref = 0.5*( Tg + Tair );	/* evaluate properties at the average temperature */
	h = h_sphere_in_air(D_GLOBE, Tref, Pair, speed);	
	v = 0.5*Qr
		- h/(STEFANB*EMIS_GLOBE)*(Tg - Tair)		
		+ term3;
	
	return v;
}

/* ============================================================================
 *  Purpose: to calculate the first derivative of the right-hand side of Eq.(17).
 *	
 *  Author:  NIES
 */

static double d_Tglobe_rhs(double Tg, double Tair, double Pair,
						   double speed, double Qr, double term3)
{
	double v;

	v = Tglobe_rhs(Tg + DELTA_DERIV, Tair, Pair, speed, Qr, term3);
	v -= Tglobe_rhs(Tg - DELTA_DERIV, Tair, Pair, speed, Qr, term3);
	v /= (2 * DELTA_DERIV);
	
	return v;
}
