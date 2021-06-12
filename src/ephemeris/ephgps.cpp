#include "navigation/ephemeris/ephgps.h"
#include "navigation/timemodule.h"
#include "navigation/matrix.h"

CEphGps::CEphGps() {
    etype_.SetEllipsoidParam(WGS84);
}

bool CEphGps::satpos(Sattime time, sat_s &sat) {
    if(sat._eph->_sys != sat._obs->_sys || sat._eph->_prn != sat._obs->_sat) return false;
    // refer to GPS ICD
    Sattime t_ = time - sat._eph->_toe;
    double tk = t_._Week * 604800.0 + t_._Sow;
    // if (!isExpired(tk)) return false;
    double n0 = Calculaten0(*sat._eph);
    double n = n0 + sat._eph->_Deltan, Mk = sat._eph->_M0 + n * tk;
    double Ek = CalculateEk(Mk, *sat._eph), vk = CalculateVk(Ek, *sat._eph);
    double Phik = sat._eph->_Omega + vk;
    double Phik2 = Phik * 2.0, sinPhik2 = sin(Phik2), cosPhik2 = cos(Phik2);
    double deltauk = sat._eph->_Cus * sinPhik2 + sat._eph->_Cuc * cosPhik2;
    double deltark = sat._eph->_Crs * sinPhik2 + sat._eph->_Crc * cosPhik2;
    double deltaik = sat._eph->_Cis * sinPhik2 + sat._eph->_Cic * cosPhik2;
    double uk = Phik + deltauk;
    double rk = sat._eph->_sqrtA * sat._eph->_sqrtA * (1 - sat._eph->_ecc * cos(Ek)) + deltark;
    double ik = sat._eph->_I0 + deltaik + sat._eph->_Idot * tk;
    double orbpos[2];
    CalOrbPos(rk, uk, orbpos);
    double omegak = sat._eph->_Omega0 + (sat._eph->_Omega_dot - etype_._rotation) *
        tk - etype_._rotation * sat._eph->_toe._Sow;
    CalEcefPos(orbpos, omegak, ik, (sat._pos));
    return true;
}

double CEphGps::Calculaten0(nav_t nav) {
    double A = nav._sqrtA * nav._sqrtA;
    return sqrt(etype_._miu / (A * A * A));
}

double CEphGps::CalculateEk(double Mk, nav_t nav) {
    double Ek0 = Mk, Ek = 0;
    int count = 0;
    while (count <= 10) {
        Ek = Mk + nav._ecc * sin(Ek0);
        double error = abs(Ek - Ek0);
        if(error < 1e-13)
            break;
        Ek0 = Ek; count ++;
    }
    return Ek;
}

double CEphGps::CalculateVk(double Ek, nav_t nav) {
    double e = nav._ecc, e2 = nav._ecc * nav._ecc;
    double sinEk = sin(Ek), cosEk = cos(Ek);
    double up = sqrt(1 - e2) * sinEk;
    double down = cosEk - e;
    return atan2(up, down);
}

void CEphGps::CalOrbPos(double rk, double uk, double* orbpos) {
    orbpos[0] = rk * cos(uk);
    orbpos[1] = rk * sin(uk);
}

void CEphGps::CalEcefPos(double* orbpos, double omegak, double ik, double* ecefpos) {
    double cosOmegak = cos(omegak), sinOmegak = sin(omegak);
    double cosIk = cos(ik), sinIk = sin(ik);
    ecefpos[0] = orbpos[0] * cosOmegak - orbpos[1] * cosIk * sinOmegak;
    ecefpos[1] = orbpos[0] * sinOmegak + orbpos[1] * cosIk * cosOmegak;
    ecefpos[2] = orbpos[1] * sinIk;
}

bool CEphGps::isExpired(double tk) {
    if (abs(tk) <= 2* 3600)
        return true;
    return false;
}

bool CEphGps::relfix(Sattime time, sat_s &sat) {
    Sattime t_ = time - sat._eph->_toe;
    double tk = t_._Week * 604800.0 + t_._Sow;
    if (!isExpired(tk)) return false;
    double n0 = Calculaten0(*sat._eph);
    double n = n0 + sat._eph->_Deltan, Mk = sat._eph->_M0 + n * tk;
    double Ek = CalculateEk(Mk, *sat._eph), vk = CalculateVk(Ek, *sat._eph);
    double F = -2.0 * sqrt(etype_._miu) / (VEL_LIGHT * VEL_LIGHT);
    double rel = F * sat._eph->_ecc * sat._eph->_sqrtA * sin(Ek);
    sat._clk[0] += rel;

    double Ekdot = n / (1.0 - sat._eph->_ecc * cos(Ek));
    rel = F * sat._eph->_ecc * sat._eph->_sqrtA * cos(Ek) * Ekdot;
    return true;
}

void CEphGps::earthRotateFix(double deltat, sat_s &sat) {
    double angle = deltat * etype_._rotation;
    Matrix<double, 3, 1> pos, result;
    Matrix<double, 3, 3> R;
    for (int i = 0; i < 3; ++i)
        pos(i, 0) = sat._pos[i]; 

    R(0, 0) = cos(angle);  R(0, 1) = sin(angle); R(0, 2) = 0;
    R(1, 0) = -sin(angle); R(1, 1) = cos(angle); R(1, 2) = 0;
    R(2, 0) = 0;           R(2, 1) = 0;          R(2, 2) = 1;
    result = R * pos;

    for (int i = 0; i < 3; ++i)
        sat._pos[i] = result(i, 0); 
}