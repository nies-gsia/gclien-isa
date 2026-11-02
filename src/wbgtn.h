/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Max Lieblich
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
 * Derived from WBGT, Version 1.1, Copyright (c) 2008, UChicago Argonne, LLC
 * (James C. Liljegren). See wbgtn.cpp and the LICENSE file for the original
 * license and the required acknowledgment.
 *
 * Modifications copyright (c) 2026 Jun'ya Takakura,
 * National Institute for Environmental Studies (NIES).
 * The modifications made by NIES are licensed under the BSD 3-Clause License
 * (SPDX-License-Identifier: BSD-3-Clause; see the LICENSE file).
 * NIES modifications: #define constants and max/min macros replaced with
 * constexpr and std::max/std::min; declarations updated for C++.
 */

#ifndef WBGT_H
#define WBGT_H

/*
*  mathematical constants
*/
constexpr double PI            = 3.1415926535897932;
constexpr double TWOPI         = 6.2831853071795864;
constexpr double DEG_RAD       = 0.017453292519943295;
constexpr double RAD_DEG       = 57.295779513082323;
/*
*  physical constants
*/
constexpr double SOLAR_CONST   = 1367.;
constexpr double GRAVITY       = 9.807;
constexpr double STEFANB       = 5.6696E-8;
constexpr double Cp            = 1003.5;
constexpr double M_AIR         = 28.97;
constexpr double M_H2O         = 18.015;
constexpr double RATIO         = Cp * M_AIR / M_H2O;
constexpr double R_GAS         = 8314.34;
constexpr double R_AIR         = R_GAS / M_AIR;
constexpr double Pr            = Cp / (Cp + 1.25 * R_AIR);
/*
*  wick constants
*/
constexpr double EMIS_WICK     = 0.95;
constexpr double ALB_WICK      = 0.4;
constexpr double D_WICK        = 0.007;
constexpr double L_WICK        = 0.0254;
/*
*  globe constants
*/
constexpr double EMIS_GLOBE    = 0.95;
constexpr double ALB_GLOBE     = 0.05;
constexpr double D_GLOBE       = 0.0508;
/*
*  surface constants
*/
constexpr double EMIS_SFC      = 0.999;
constexpr double ALB_SFC       = 0.45;
/*
*  computational and physical limits
*/
constexpr double CZA_MIN       = 0.00873;
constexpr double NORMSOLAR_MAX = 0.85;
constexpr double REF_HEIGHT    = 2.0;
constexpr double MIN_SPEED     = 0.13;
constexpr double CONVERGENCE   = 0.02;
constexpr int    MAX_ITER      = 500;
/*
*  computational parameter
*/
constexpr double DELTA_DERIV   = 1e-10;

/*
 * All function declarations in this file are extracted from Liljegren's original code
 * in wbgt.c and converted to ANSI declaration style from K&R style.
 * Old K&R variable declarations are copied here to document the code.
 */

/*
int calc_wbgt(int year, int hour, int day, int month, int minute, int gmt,
              int avg, double lat, double lon, double solar, double pres, double Tair,
              double relhum, double speed, double zspeed, double dT, int urban,
              double* est_speed, double* Tg, double* Tnwb, double* Tpsy, double* Twbg);
*/

int calc_wbgt(int year, int month, int day, int hour, int minute, int gmt,
    int avg, double lat, double lon, double solar, double pres, double Tair,
    double relhum, double speed, double zspeed, double dT, int urban,
    double* est_speed, double* Tg, double* Tnwb, double* Tpsy, double* Twbg);


/*
 int	year,	4-digit, e.g. 2007
month,	month (1-12) or month = 0 implies iday is day of year
day,		day of month or day of year (1-366)
hour,		hour in local standard time (LST)
minute,	minutes past the hour
gmt,		LST-GMT difference, hours (negative in USA)
avg,		averaging time of meteorological inputs, minutes
urban;	select "urban" (1) or "rural" (0) wind speed power law exponent

double	lat,	north latitude, decimal
lon,		  east longitude, decimal (negative in USA)
solar,	  solar irradiance, W/m2
pres,		  barometric pressure, mb
Tair,		  air (dry bulb) temperature, degC
relhum,	  relative humidity, %
speed,	  wind speed, m/s
zspeed,	  height of wind speed measurement, m
dT,		    vertical temperature difference (upper minus lower), degC

*est_speed,	 estimated speed at reference height, m/s
*Tg,		  globe temperature, degC
*Tnwb,	  natural wet bulb temperature, degC
*Tpsy,	  psychrometric wet bulb temperature, degC
*Twbg;	  wet bulb globe temperature, degC
*/


int	calc_solar_parameters(int year, int month, double day, double lat, double lon,
                          double *solar, double *cza, double *fdir);
/*
 int	year,		 4-digit year, e.g., 2007
month;	 2-digit month; month = 0 implies day = day of year

double day;	day.fraction of month if month > 0;
else day.fraction of year if month = 0 (GMT)
double	lat,	 north latitude
lon,		   east latitude (negative in USA)
*solar,	   solar irradiance (W/m2)
*cza,		   cosine of solar zenith angle
*fdir;	   fraction of solar irradiance due to direct beam
*/

int	calc_solar_cza(int year, int month, double day, double lat, double lon, double* cza);

double Twb(double Tair, double rh, double Pair, double speed,
          double solar, double fdir, double cza, int rad);
/*
double Tair,		 air (dry bulb) temperature, degC
rh,
Pair,		barometric pressure, mb
speed,	wind speed, m/s
solar,	solar irradiance, W/m2
fdir,		fraction of solar irradiance due to direct beam
cza;		cosine of solar zenith angle

int	rad;		switch to enable/disable radiative heating;
            no radiative heating --> pyschrometric wet bulb temp
*/

double h_cylinder_in_air(double diameter, double length, double Tair, double Pair, double speed);
/*
  double	diameter,	cylinder diameter, m
length,	 cylinder length, m
Tair,		 air temperature, K
Pair,		 barometric pressure, mb
speed;	 fluid (wind) speed, m/s
*/

double Tglobe(double Tair, double rh, double Pair, double speed,
             double solar, double fdir, double cza);
/*
  double Tair,		 air (dry bulb) temperature, degC
rh,		 relative humidity, fraction between 0 and 1
Pair,		 barometric pressure, mb
speed,	 wind speed, m/s
solar,	 solar irradiance, W/m2
fdir,		 fraction of solar irradiance due to direct beam
cza;		 cosine of solar zenith angle
*/

double h_sphere_in_air(double diameter, double Tair, double Pair, double speed);
  /*
  double	diameter,	 sphere diameter, m
  Tair,		 air temperature, K
  Pair,		 barometric pressure, mb
  speed;	 fluid (air) speed, m/s
  */

double esat(double tk, int phase);
/*
double	tk;	air temperature, K
int	phase;
*/

double dew_point(double e, int phase);
/*
double	e;	vapor pressure, mb
int	phase;
*/

double viscosity(double Tair);
/* double	Tair;	air temperature, K */

double thermal_cond(double Tair);
/*double	Tair;	air temperature, K */

double diffusivity(double Tair, double Pair);
/*
double	Tair,	Air temperature, K
Pair; Barometric pressure, mb
*/

double evap(double Tair);
/* double	Tair;	air temperature, K */

double emis_atm(double Tair, double rh);
/*
double	Tair,	air temperature, K
rh;	relative humidity, fraction between 0 and 1
*/

int solarposition(int year, int month, double day, double days_1900, double latitude,
                  double longitude, double *ap_ra, double *ap_dec, double *altitude,
                  double *refraction, double *azimuth, double *distance);
/*
int year, Four digit year (Gregorian calendar).
          [1950 through 2049; 0 o.k. if using days_1900]
month;    Month number.
          [1 through 12; 0 o.k. if using daynumber for day]
double day, Calendar day.fraction, or daynumber.fraction.
            [If month is NOT 0: 0 through 32; 31st @ 18:10:00 UT = 31.75694
             If month IS 0: 0 through 367; 366 @ 18:10:00 UT = 366.75694]
days_1900, Days since 1900 January 0 @ 00:00:00 UT.
    [18262.0 (1950/01/00) through 54788.0 (2049/12/32);
     1990/01/01 @ 18:10:00 UT = 32873.75694;
     0.0 o.k. if using {year, month, day} or
     {year, daynumber}]
latitude, Observation site geographic latitude.
    [degrees.fraction, North positive]
longitude,     Observation site geographic longitude.
    [degrees.fraction, East positive]
*ap_ra,        Apparent solar right ascension.
    [hours; 0.0 <= *ap_ra < 24.0]
*ap_dec,       Apparent solar declination.
    [degrees; -90.0 <= *ap_dec <= 90.0]
*altitude,     Solar altitude, uncorrected for refraction.
    [degrees; -90.0 <= *altitude <= 90.0]
*refraction,   Refraction correction for solar altitude.
  Add this to altitude to compensate for refraction.
    [degrees; 0.0 <= *refraction]
*azimuth, Solar azimuth.
    [degrees; 0.0 <= *azimuth < 360.0, East is 90.0]
*distance; Distance of Sun from Earth (heliocentric-geocentric).
    [astronomical units; 1 a.u. is mean distance]
 */

int daynum(int year, int month, int day);
/*  int year, month, day;*/

double est_wind_speed(double speed, double zspeed, int stability_class, int urban);
/*
int	stability_class, urban;
double	speed, zspeed;
*/

int stab_srdt(int daytime, double speed, double solar, double dT);
/*
int	daytime;
double	speed,
solar,
dT;
*/


#endif /* WBGT_H */
