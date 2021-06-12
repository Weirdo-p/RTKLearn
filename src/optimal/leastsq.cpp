#include "navigation/optimal/leastsq.h"
#include "navigation/pnt/pntbase.h"
#include "navigation/pnt/pntrtk.h"

bool CLeastsq::optimize(MatrixXd B, MatrixXd P, MatrixXd w, int nobs) {
    int istrue = 0;
    MatrixXd inv = B.transpose() * P * B;
    inv = inv.inverse(istrue);
    _Q = inv;
    _x = inv * B.transpose() * P * w;

    MatrixXd v = B * _x ;
    v = v - w;
    auto sigma = v.transpose() * P * v;
    _sigma = sigma(0, 0) / (nobs - _opt->_nsys - 3);
    return bool(istrue);
}

bool CLeastsq::optimizeAVD(MatrixXd B, MatrixXd P, MatrixXd w, int nobs) {
    int istrue = 0;
    MatrixXd inv = B.transpose() * P * B;
    inv = inv.inverse(istrue);
    _vel_Q = inv;
    _vel = inv * B.transpose() * P * w;
    MatrixXd v = B * _x ;
    v = v - w;
    auto sigma = v.transpose() * P * v;
    _vel_sigma = sigma(0, 0) / (nobs - _opt->_nsys - 3);
    return bool(istrue);
}


MatrixXd CLeastsq::Getx() {
    return _x;
}

CLeastsq::CLeastsq() {

}

CLeastsq::CLeastsq(prcopt* opt) {
    if (!_opt) {
        _opt = new prcopt;
        memcpy(_opt, opt, sizeof(prcopt));
    }
}

void CLeastsq::getl(sat* sats, double* basepos, int* refsats, MatrixXd pos, MatrixXd &w) {
    int nobs = sats->_nsats;
    int ref_prn = -1, sysflag = -1;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
    #pragma omp parallel for
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1]._sat[isat]._isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((_opt->_navsys & SYS_ARRAY[i]) == sats[1]._sat[isat]._sys) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1]._sat[isat]._prn == ref_prn && sats[1]._sat[isat]._sys == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        double dist_base_ref = CPntbase::getSatUserPos(ref_sat_b, basepos);
        double dist_base_com = CPntbase::getSatUserPos(sats[0]._sat[isat], basepos);
        double dist_rover_ref = CPntbase::getSatUserPos(ref_sat_r, rover_pos);
        double dist_rover_com = CPntbase::getSatUserPos(sats[1]._sat[isat], rover_pos);
        double dist_cal = dist_rover_com - dist_rover_ref - dist_base_com + dist_base_ref;
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sats[1]._sat[isat]._obs->_L[i] - ref_sat_r._obs->_L[i]-
                    sats[0]._sat[isat]._obs->_L[i] + ref_sat_b._obs->_L[i];
                double pseu_obs = sats[1]._sat[isat]._obs->_P[i] - ref_sat_r._obs->_P[i]- 
                    sats[0]._sat[isat]._obs->_P[i] + ref_sat_b._obs->_P[i];
                double freq = CPntbase::GetFreq(sysflag, FREQ_ARRAY[i]);
                w(num ++, 0) = (pseu_obs - dist_cal);
                w(num ++, 0) = (phase_obs * (VEL_LIGHT / freq) - dist_cal);
            }
        }
    }
}

bool CLeastsq::optimize(sat* sats_epoch, res_t &res) {
    int row = 0, col = 0, iter = 0;
    const int MAXITER = 10;
    int refsats[MAXSYS];
    int nobs = CPntrtk::obsnumber(sats_epoch);
    chooseref(sats_epoch[1], refsats);

    /* LEAST SQUARE METHOD */
    MatrixXd pos;
    getDesignDim(sats_epoch[1], nobs, row, col);
    while (iter < MAXITER) {
        MatrixXd B(row, col);
        MatrixXd w(row, 1);
        MatrixXd P(row, row);
        if (iter == 0) {
            pos.resize(col, 1); pos.Zero();
            for (int i = 0; i < 3; ++i) pos(i, 0) = res._rpos_ecef[i];
        }
        double rover_pos[3] = {0};
        
        for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
        B.Zero(); P.Zero(); w.Zero();
        getDesign(sats_epoch, nobs, rover_pos, refsats, B);
        getl(sats_epoch, _opt->_base, refsats, pos, w);
        getweight(sats_epoch, refsats, nobs, P);
        int flag = 0;
        if(!optimize(B, P, w, nobs)) {
            continue;
        }

        MatrixXd x_pos(3, 1);
        for (int i = 0; i < _x.row(); ++i) {
            if (i < 3) {
                pos(i, 0) += _x(i, 0); x_pos(i, 0) = _x(i, 0);
            } else 
                pos(i, 0) = _x(i, 0);
        }
        _x = pos;
        if (x_pos.norm() < 1E-6) break;
        iter ++;
    }
    for (int i = 0; i < 3; ++ i) res._rpos_ecef[i] = pos(i, 0);
}

MatrixXd CLeastsq::GetVar() {
    return _Q;
}

double CLeastsq::GetInternalSigma() {
    return _sigma;
}

MatrixXd CLeastsq::GetState() {
    return _x;
}

MatrixXd CLeastsq::GetVel() {
    return _vel;
}

MatrixXd CLeastsq::GetVelVar() {
    return _vel_Q;
}

double CLeastsq::GetVelInternalSigma() {
    return _vel_sigma;
}