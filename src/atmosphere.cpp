#include "navigation/atmosphere.h"
#include "navigation/coors.h"
#include <math.h>

double Hopefield(const double E, const double H, double t0,
                 double p0, double RH0) {
    double E1 = Rad2Deg(E);
    double hw = 11000, H0 = 0.0;
    double T0 = t0 + 273.16;
    double RH = RH0 * exp(-0.0006396 * (H - H0));
    double p = p0 * pow(1 - 0.0000226 * (H - H0), 5.225);
    double T = T0 - 0.0065 * (H - H0);
    double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
    double hd = 40136 + 148.72 * (T0 - 273.16);
    double Kd = 155.2e-7 * p / T * (hd - H);
    double Kw = 155.2e-7 * 4810.0 / (T * T) * e * (hw - H);
    double left = Kd / sin(Deg2Rad(sqrt(E1 * E1 + 6.25)));
    double right = Kw / sin(Deg2Rad(sqrt(E1 * E1 + 2.25)));
    return left + right;
}

double CAtmosphere::klobuchar(double* pos, sat_s sat) {

    const double ion[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double CLIGHT = VEL_LIGHT;

    Sattime gpst = sat._obs->_time;

    double azel[2] = {sat._azi, sat._elev};
    double tt, f, psi, phi, lam, amp, per, x;
    int week;

    if (pos[2] < -1E3 || azel[1] <= 0)
    {
        return 0.0;
    }


    /* earth centered angle (semi-circle) */
    psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;

    /* subionospheric latitude/longitude (semi-circle) */
    phi = pos[0] / PI + psi * cos(azel[0]);
    if (phi > 0.416) phi = 0.416;
    else if (phi < -0.416) phi = -0.416;
    lam = pos[1] / PI + psi * sin(azel[0]) / cos(phi * PI);

    /* geomagnetic latitude (semi-circle) */
    phi += 0.064 * cos((lam - 1.617) * PI);

    /* local time (s) */
    //tt = 43200.0 * lam + time2gpst(t, &week);
    tt = 43200.0 * lam + gpst._Sow;
    //std::cout << "Sec of week:" << gpst.sec_of_week;
    week = gpst._Week;

    tt -= floor(tt / 86400.0) * 86400.0; /* 0<=tt<86400 */

    /* slant factor */
    f = 1.0 + 16.0 * pow(0.53 - azel[1] / PI, 3.0);

    /* ionospheric delay */
    amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]));
    per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]));
    amp = amp < 0.0 ? 0.0 : amp;
    per = per < 72000.0 ? 72000.0 : per;
    x = 2.0 * PI * (tt - 50400.0) / per;

    return CLIGHT * f * (fabs(x) < 1.57 ? 5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0)) : 5E-9);
}
