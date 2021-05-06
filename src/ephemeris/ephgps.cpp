#include "navigation/ephemeris/ephgps.h"
#include "navigation/timemodule.h"
#include "matrix.h"

CEphGps::CEphGps() {
    etype_.SetEllipsoidParam(WGS84);
}

bool CEphGps::satpos(Sattime time, sat_s &sat) {
    if(sat.eph_->sys_ != sat.obs_->sys || sat.eph_->prn_ != sat.obs_->sat) return false;
    // refer to GPS ICD
    Sattime t_ = time - sat.eph_->toe_;
    double tk = t_.Week_ * 604800.0 + t_.Sow_;
    // if (!isExpired(tk)) return false;
    double n0 = Calculaten0(*sat.eph_);
    double n = n0 + sat.eph_->Deltan_, Mk = sat.eph_->M0_ + n * tk;
    double Ek = CalculateEk(Mk, *sat.eph_), vk = CalculateVk(Ek, *sat.eph_);
    double Phik = sat.eph_->Omega_ + vk;
    double Phik2 = Phik * 2.0, sinPhik2 = sin(Phik2), cosPhik2 = cos(Phik2);
    double deltauk = sat.eph_->Cus_ * sinPhik2 + sat.eph_->Cuc_ * cosPhik2;
    double deltark = sat.eph_->Crs_ * sinPhik2 + sat.eph_->Crc_ * cosPhik2;
    double deltaik = sat.eph_->Cis_ * sinPhik2 + sat.eph_->Cic_ * cosPhik2;
    double uk = Phik + deltauk;
    double rk = sat.eph_->sqrtA_ * sat.eph_->sqrtA_ * (1 - sat.eph_->ecc_ * cos(Ek)) + deltark;
    double ik = sat.eph_->I0_ + deltaik + sat.eph_->Idot_ * tk;
    double orbpos[2];
    CalOrbPos(rk, uk, orbpos);
    double omegak = sat.eph_->Omega0_ + (sat.eph_->Omega_dot_ - etype_.rotation_) *
        tk - etype_.rotation_ * sat.eph_->toe_.Sow_;
    CalEcefPos(orbpos, omegak, ik, (sat.pos_));
    return true;
}

double CEphGps::Calculaten0(nav_t nav) {
    double A = nav.sqrtA_ * nav.sqrtA_;
    return sqrt(etype_.miu_ / (A * A * A));
}

double CEphGps::CalculateEk(double Mk, nav_t nav) {
    double Ek0 = Mk, Ek = 0;
    int count = 0;
    while (count <= 10) {
        Ek = Mk + nav.ecc_ * sin(Ek0);
        double error = abs(Ek - Ek0);
        if(error < 1e-13)
            break;
        Ek0 = Ek; count ++;
    }
    return Ek;
}

double CEphGps::CalculateVk(double Ek, nav_t nav) {
    double e = nav.ecc_, e2 = nav.ecc_ * nav.ecc_;
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
    Sattime t_ = time - sat.eph_->toe_;
    double tk = t_.Week_ * 604800.0 + t_.Sow_;
    if (!isExpired(tk)) return false;
    double n0 = Calculaten0(*sat.eph_);
    double n = n0 + sat.eph_->Deltan_, Mk = sat.eph_->M0_ + n * tk;
    double Ek = CalculateEk(Mk, *sat.eph_), vk = CalculateVk(Ek, *sat.eph_);
    double F = -2.0 * sqrt(etype_.miu_) / (VEL_LIGHT * VEL_LIGHT);
    double rel = F * sat.eph_->ecc_ * sat.eph_->sqrtA_ * sin(Ek);
    sat.clk[0] += rel;

    double Ekdot = n / (1.0 - sat.eph_->ecc_ * cos(Ek));
    rel = F * sat.eph_->ecc_ * sat.eph_->sqrtA_ * cos(Ek) * Ekdot;
    return true;
}

void CEphGps::earthRotateFix(double deltat, sat_s &sat) {
    double angle = deltat * etype_.rotation_;
    Matrix<double, 3, 1> pos, result;
    Matrix<double, 3, 3> R;
    for (int i = 0; i < 3; ++i)
        pos(i, 0) = sat.pos_[i]; 

    R(0, 0) = cos(angle);  R(0, 1) = sin(angle); R(0, 2) = 0;
    R(1, 0) = -sin(angle); R(1, 1) = cos(angle); R(1, 2) = 0;
    R(2, 0) = 0;           R(2, 1) = 0;          R(2, 2) = 1;
    result = R * pos;

    for (int i = 0; i < 3; ++i)
        sat.pos_[i] = result(i, 0); 
}