#include "navigation/pnt/pntbase.h"
#include "navigation/rinex/rinex304.h"
#include "navigation/ephemeris/ephgps.h"
#include "navigation/ephemeris/ephbds.h"
#include "navigation/coors.h"
#include "navigation/atmosphere.h"
#include "navigation/ambiguity/lambda.h"
#include "navigation/timemodule.h"
#include <stdio.h>
#include <string.h>
#include <omp.h>

CPntbase::CPntbase() {
    _rnx = nullptr; _satpos = nullptr;
    _opt = nullptr; _optimizer = nullptr;
    _res = new res_t;
    memset(_res, 0, sizeof(res_t));
}

CPntbase::CPntbase(prcopt opt) {
    _opt = new prcopt;
    memcpy(_opt, &opt, sizeof(opt));
    _rnx = nullptr; _satpos = nullptr;
    _res = new res_t;
    memset(_res, 0, sizeof(res_t));

    if (_optimizer)
        return;
    if (opt._proctype == PROC_LS)
        _optimizer = new CLeastsq();
    else if (opt._proctype == PROC_KF)
        _optimizer = new CRtkekf();
}

void CPntbase::readRinex(char** paths, int n) {
    bindRinex(paths[0]);
    for (int i = 0; i < n; ++ i) {
        _rnx->setsites(_opt->_sitenum);
        _rnx->decode(paths[i], *_opt);
    }
}

void CPntbase::bindRinex(char* path) {
    string ver = this->_rnx->readver(path);
    if (ver == "3.04")
        _rnx = new CDecodeRnx304;
}

CPntbase::~CPntbase() {
    if (_rnx) {
        delete _rnx; _rnx = nullptr;
    }
    if (_res) {
        delete _res; _res = nullptr;
    }
    if (_satpos){
        delete _satpos; _satpos = nullptr;
    }
    if(_opt) {
        delete _opt; _opt = nullptr;
    }
}

bool CPntbase::satclk(sat* sats) {
    int sitenum = _opt->_sitenum;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site]._nsats; ++ i_sat) {
            if (!sats[i_site]._sat[i_sat]._obs || !sats[i_site]._sat[i_sat]._eph)
                continue;
            bindEph(sats[i_site]._sat[i_sat]._sys);
            Sattime tobs = sats[i_site]._sat[i_sat]._obs->_time;
            double P = searchpseu(sats[i_site]._sat[i_sat]._obs->_P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site]._sat[i_sat]._sys << "prn: " << 
                        sats[i_site]._sat[i_sat]._prn << endl;
            tobs = tobs - P / VEL_LIGHT;
            memset(sats[i_site]._sat[i_sat]._clk, 0, sizeof(double) * 2);
            _satpos->satclk(tobs, sats[i_site]._sat[i_sat]);

            delete _satpos; _satpos = nullptr;
        }
    }
    return true;
}

void CPntbase::bindEph(int sys) {
    switch (sys) {
    case SYS_GPS:
        _satpos = new CEphGps;
        break;
    case SYS_BDS:
        _satpos = new CEphBds;
        break;
    default:
        _satpos = new CEphBase;
        break;
    }
}

bool CPntbase::satpos(sat* sats) {
    int sitenum = _opt->_sitenum;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site]._nsats; ++ i_sat) {
            if (!sats[i_site]._sat[i_sat]._obs || !sats[i_site]._sat[i_sat]._eph)
                continue;
            bindEph(sats[i_site]._sat[i_sat]._sys);
            Sattime tobs = sats[i_site]._sat[i_sat]._obs->_time;
            double P = searchpseu(sats[i_site]._sat[i_sat]._obs->_P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site]._sat[i_sat]._sys << "prn: " << 
                        sats[i_site]._sat[i_sat]._prn << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site]._sat[i_sat]._clk[0];
            sats[i_site]._sat[i_sat]._eph->_sig = tobs;
            _satpos->satpos(tobs, sats[i_site]._sat[i_sat]);
            delete _satpos; _satpos = nullptr;
        }
    }
    return true;
}

double CPntbase::searchpseu(double* P) {
    double p = 0;
    for (int i = 0; i < MAXFREQ; ++ i) {
        if (P[i] != 0) {
            p = P[i]; break;
        }
    }
    return p;
}

nav_t* CPntbase::searchnav(Sattime time, int sys, int prn, nav* navs) {
    int delay = 0;
    bool isfind = false;
    if (sys == SYS_GPS) delay = 7200;
    else if (sys == SYS_BDS) {
        time = time - BDT2GPST;
        delay = 3600;
    }
    // nav_t* sat_nav = nullptr;
    #pragma omp parallel for
    for (int i_nav = navs->_num; i_nav > 0; -- i_nav) {
        if (sys != navs->_msg[i_nav]._sys || prn != navs->_msg[i_nav]._prn)
            continue;
        if (abs(Sattimediff(time, navs->_msg[i_nav]._toe)) > delay)
            continue;
        isfind = true;
        return &navs->_msg[i_nav];
    }
    if (!isfind) {
        cout << "no ephemeris found: sys: " << sys << " prn: " << prn << endl;
        return nullptr;
    }
}

bool CPntbase::satvel(sat* sats) {
    int sitenum = _opt->_sitenum;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site]._nsats; ++ i_sat) {
            bindEph(sats[i_site]._sat[i_sat]._sys);
            Sattime tobs = sats[i_site]._sat[i_sat]._obs->_time;
            double P = searchpseu(sats[i_site]._sat[i_sat]._obs->_P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site]._sat[i_sat]._sys << "prn: " << 
                        sats[i_site]._sat[i_sat]._prn << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site]._sat[i_sat]._clk[0];
            _satpos->satvel(tobs, sats[i_site]._sat[i_sat]);
            delete _satpos; _satpos = nullptr;
        }
    }
    return true;
}

bool CPntbase::relativeeffect(sat* sats) {
    int sitenum = _opt->_sitenum;
    for (int i_site = 0; i_site < sitenum; ++ i_site)
        for (int i_sat = 0; i_sat < sats[i_site]._nsats; ++ i_sat) {
            bindEph(sats[i_site]._sat[i_sat]._sys);
            Sattime tobs = sats[i_site]._sat[i_sat]._obs->_time;
            double P = searchpseu(sats[i_site]._sat[i_sat]._obs->_P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site]._sat[i_sat]._sys << "prn: " << 
                        sats[i_site]._sat[i_sat]._prn << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site]._sat[i_sat]._clk[0];
            _satpos->relfix(tobs, sats[i_site]._sat[i_sat]);
            delete _satpos; _satpos = nullptr;
        }
}

void CPntbase::earthRotateFix(sat* sats) {
    int sitenum = _opt->_sitenum;
    for (int i_site = 0; i_site < sitenum; ++ i_site)
        for (int i_sat = 0; i_sat < sats[i_site]._nsats; ++ i_sat) {
            if (!sats[i_site]._sat[i_sat]._eph) continue;
            bindEph(sats[i_site]._sat[i_sat]._sys);
            Sattime dt = sats[i_site]._sat[i_sat]._obs->_time -sats[i_site]._sat[i_sat]._eph->_sig;
            _satpos->earthRotateFix(dt._2sec(), sats[i_site]._sat[i_sat]);
            delete _satpos; _satpos  = nullptr;
        }
}

void CPntbase::satazel(double* site, sat &sat) {
    int satnum = sat._nsats;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (*site == 0) {
            sat._sat[i_sat]._elev = sat._sat[i_sat]._azi = 0.0;
            break;
        } else {
            double neu[3];
            bindEph(sat._sat[i_sat]._sys);
            XYZ2NEU(site, sat._sat[i_sat]._pos, _satpos->GetElli(), neu);
            sat._sat[i_sat]._elev = atan(neu[2] / sqrt(neu[0] * neu[0] + neu[1] * neu[1]));
            sat._sat[i_sat]._azi = atan2(neu[0], neu[1]);
            delete _satpos; _satpos = nullptr;
        }
    }
}

int CPntbase::usesys(sat &sats) {
    int satnum = sats._nsats, num = 0;
    int used = 0x00;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat)
        used |= sats._sat[i_sat]._sys;
    for (int i = 0; i < MAXSYS; ++i) 
        if ((used & SYS_ARRAY[i]) == SYS_ARRAY[i])
            num ++;
    return num;
}

void CPntbase::GetLC(sat sat, MatrixXd &w) {
    int satnum = sat._nsats, num = 0;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sat._sat[i_sat]._isused) continue;
        if (!sat._sat[i_sat]._isused) continue;
        int p_pos[2] = {0};  // used pseudorange position
        for (int i = 0, j = 0; i < MAXFREQ; ++ i) 
            if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i])
                p_pos[j++] = i;
        double L1 = GetFreq(sat._sat[i_sat]._sys, FREQ_ARRAY[p_pos[0]]);
        double L2 = GetFreq(sat._sat[i_sat]._sys, FREQ_ARRAY[p_pos[1]]);
        double L1_square = L1 * L1, L2_square = L2 * L2;
        double IF = L1_square / (L1_square - L2_square) * sat._sat[i_sat]._obs->_P[p_pos[0]] -
                    L2_square / (L1_square - L2_square) * sat._sat[i_sat]._obs->_P[p_pos[1]];
        if(sat._sat[i_sat]._sys == SYS_BDS)
            IF -= VEL_LIGHT * L1_square * sat._sat[i_sat]._eph->_Tgd[0] / 
                  (L1_square - L2_square);
        w(num , 0) += IF; num ++;
    }
}

void CPntbase::GetNonCombine(sat sat, MatrixXd &w) {
    int satnum = sat._nsats, num = 0;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sat._sat[i_sat]._isused) continue;
        int p_pos;  // used pseudorange position
        for (int i = 0; i < MAXFREQ; ++ i) 
            if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i])
                p_pos = i;
        double L1 = sat._sat[i_sat]._obs->_P[p_pos];
        w(num, 0) += L1; num ++;
    }
}


void CPntbase::setopt(prcopt opt) {
    // memcpy(&_opt, &opt, sizeof(opt));
    if (_opt) {
        _opt = &opt;
    } else {
        memcpy(_opt, &opt, sizeof(opt));
    }
}

double CPntbase::GetFreq(int sys, int freqflag) {
    switch (sys) {
    case SYS_GPS:
        if (freqflag == FREQTYPE_L1)
            return FREQ_L1;
        else if (freqflag == FREQTYPE_L2)
            return FREQ_L2;
        else return FREQ_L5;
    case SYS_BDS:
        if (freqflag == FREQTYPE_L1)
            return FREQ_B1I;
        else if (freqflag == FREQTYPE_L2)
            return FREQ_B3I;
        else return FREQ_B2;
    default:
        return 0;
    }
}

void CPntbase::tropfix(sat& sats, MatrixXd &w, double H) {
    int satnum = sats._nsats, num = 0;
    for (int isat = 0; isat < satnum; ++ isat) {
        if (!sats._sat[isat]._isused) continue;
        double T = Hopefield(sats._sat[isat]._elev, H);
        w(num, 0) -= T; num ++;
    }
}

void CPntbase::setres(double* sitepos_ecef, double* sitepos_blh, Ellipsoid type, MatrixXd result) {
    for (int i = 0; i < 3; ++ i) {
        sitepos_ecef[i] = result(i, 0);
    }
    XYZ2BLH(sitepos_ecef, type, sitepos_blh);
}

void CPntbase::selectPos(double** xzy, double** blh, int i_site) {
    if (i_site == 0) {
        if (_res->_bpos_ecef[0] != 0) {
            *xzy = _res->_bpos_ecef;
            *blh = _res->_bpos_blh;
        } else if (_opt->_base[0] != 0) {
            memcpy(_res->_bpos_ecef, _opt->_base, sizeof(double) * 3);
            *xzy = _res->_bpos_ecef;
            *blh = _res->_bpos_blh;
        } else {
            memset(_res->_bpos_ecef, 0, sizeof(double) * 3);
            *xzy = _res->_bpos_ecef;
            *blh = _res->_bpos_blh;
        }
    } else {
        if (_res->_rpos_ecef[0] != 0) {
            *xzy = _res->_rpos_ecef;
            *blh = _res->_rpos_blh;
        } else if (_opt->_rover[0] != 0) {
            memcpy(_res->_rpos_ecef, _opt->_rover, sizeof(double) * 3);
            *xzy = _res->_rpos_ecef;
            *blh = _res->_rpos_blh;
        } else {
            memset(_res->_rpos_ecef, 0, sizeof(double) * 3);
            *xzy = _res->_rpos_ecef;
            *blh = _res->_rpos_blh;
        }
    }
}

void CPntbase::getSysObs(sat sats, int* nobs) {
    if(!nobs) return;
    int nobs_sys[MAXSYS] = {0};
    for (int isat = 0; isat < sats._nsats; ++isat) 
        for (int i = 0; i < MAXSYS; ++i) 
            if (sats._sat[isat]._sys == SYS_ARRAY[i] && sats._sat[isat]._isused)
                nobs_sys[i] += 1;
    memcpy(nobs, nobs_sys, sizeof(int) * MAXSYS);
}

double CPntbase::getSatUserPos(sat_s sats, double* sitepos) {
    double coef[3] = {0}, dist = 0;
    for (int i = 0; i < 3; ++i) {
        coef[i] = sitepos[i] - sats._pos[i];
        dist += coef[i] * coef[i];
    }
    return sqrt(dist);
}

double CPntbase::weightbyelev(double elev, double sigma0, double alpha) {
    double sigma02 = sigma0 * sigma0;
    double cosE = cos(elev);
    return (sigma02 * (1.0 * alpha * cosE * cosE));
}

void CPntbase::fixambi() {
    const int solu_num = 2;
    _res->_ambi_flag = FLOAT_SOLU;
    if (_opt->_soltype == SOLTYPE_FLOAT) return;
    auto state = _optimizer->GetState();
    auto var = _optimizer->GetVar();

    if (abs(state(0, 0)) < 1e-4)
        return;
    // double sigma = _optimizer->GetInternalSigma(); // var = var * sigma;

    int n_ambi = state.row() - 3;
    auto float_solu = state.block<3, 0>(n_ambi, 1);
    auto float_solu_var = var.block<3, 3>(n_ambi, n_ambi);
    double* fix_tmp = new double[n_ambi * solu_num];
    MatrixXd fix(n_ambi, 1);
    MatrixXd s(1, solu_num);
    lambda(n_ambi, solu_num, float_solu._2array(), float_solu_var._2array(), fix_tmp, s._2array());
    for (int i = 0; i < n_ambi; ++ i) fix(i, 0) = fix_tmp[i];
    _res->_ratio = s(0, 1) / s(0, 0);
    if (_res->_ratio < RATIO_THRES)
        return;

    _res->_ambi_flag = FIX_SOLU;
    int is_inverse = 0;
    // update result
    auto Qba = var.block<0, 3> (3, n_ambi);
    auto Qab = var.block<3, 0>(n_ambi, 3);
    auto Qbb = var.block<0, 0> (3, 3);
    auto state_float = state.block<0, 0>(3, 1);
    MatrixXd state_fix = state_float - Qba * float_solu_var.inverse(is_inverse) * (float_solu - fix);
    for (int i = 0; i < 3; ++ i) _res->_rpos_ecef[i] = state_fix(i, 0);
    delete[] fix_tmp; fix_tmp = nullptr;
}
