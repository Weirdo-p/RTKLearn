#include "navigation/pnt/pntbase.h"
#include "navigation/rinex/rinex304.h"
#include "navigation/ephemeris/ephgps.h"
#include "navigation/ephemeris/ephbds.h"
#include "navigation/coors.h"
#include "navigation/atmosphere.h"
#include "navigation/ambiguity/lambda.h"
#include <stdio.h>
#include <string.h>
#include <omp.h>

CPntbase::CPntbase() {
    rnx_ = nullptr; satpos_ = nullptr;
    opt_ = nullptr; optimizer_ = nullptr;
    res_ = new res_t;
    memset(res_, 0, sizeof(res_t));
}

CPntbase::CPntbase(prcopt opt) {
    opt_ = new prcopt;
    memcpy(opt_, &opt, sizeof(opt));
    rnx_ = nullptr; satpos_ = nullptr;
    res_ = new res_t;
    memset(res_, 0, sizeof(res_t));

    if (optimizer_)
        return;
    if (opt.proctype_ == PROC_LS)
        optimizer_ = new CLeastsq();
    else if (opt.proctype_ == PROC_KF)
        optimizer_ = new CRtkekf();
}

void CPntbase::readRinex(char** paths, int n) {
    bindRinex(paths[0]);
    for (int i = 0; i < n; ++ i) {
        rnx_->setsites(opt_->sitenum_);
        rnx_->decode(paths[i], *opt_);
    }
}

void CPntbase::bindRinex(char* path) {
    string ver = this->rnx_->readver(path);
    if (ver == "3.04")
        rnx_ = new CDecodeRnx304;
}

CPntbase::~CPntbase() {
    if (rnx_) {
        delete rnx_; rnx_ = nullptr;
    }
    if (res_) {
        delete res_; res_ = nullptr;
    }
    if (satpos_){
        delete satpos_; satpos_ = nullptr;
    }
    if(opt_) {
        delete opt_; opt_ = nullptr;
    }
}

bool CPntbase::satclk(sat* sats) {
    int sitenum = opt_->sitenum_;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site].nsats_; ++ i_sat) {
            bindEph(sats[i_site].sat_[i_sat].sys_);
            Sattime tobs = sats[i_site].sat_[i_sat].obs_->time;
            double P = searchpseu(sats[i_site].sat_[i_sat].obs_->P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site].sat_[i_sat].sys_ << "prn: " << 
                        sats[i_site].sat_[i_sat].prn_ << endl;
            tobs = tobs - P / VEL_LIGHT;
            memset(sats[i_site].sat_[i_sat].clk, 0, sizeof(double) * 2);
            satpos_->satclk(tobs, sats[i_site].sat_[i_sat]);

            delete satpos_; satpos_ = nullptr;
        }
    }
    return true;
}

void CPntbase::bindEph(int sys) {
    switch (sys) {
    case SYS_GPS:
        satpos_ = new CEphGps;
        break;
    case SYS_BDS:
        satpos_ = new CEphBds;
        break;
    default:
        satpos_ = new CEphBase;
        break;
    }
}

bool CPntbase::satpos(sat* sats) {
    int sitenum = opt_->sitenum_;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site].nsats_; ++ i_sat) {
            if (!sats[i_site].sat_[i_sat].eph_) continue;
            bindEph(sats[i_site].sat_[i_sat].sys_);
            Sattime tobs = sats[i_site].sat_[i_sat].obs_->time;
            double P = searchpseu(sats[i_site].sat_[i_sat].obs_->P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site].sat_[i_sat].sys_ << "prn: " << 
                        sats[i_site].sat_[i_sat].prn_ << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site].sat_[i_sat].clk[0];
            sats[i_site].sat_[i_sat].eph_->sig_ = tobs;
            satpos_->satpos(tobs, sats[i_site].sat_[i_sat]);
            delete satpos_; satpos_ = nullptr;
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

bool CPntbase::satvel(sat* sats) {
    int sitenum = opt_->sitenum_;
    for (int i_site = 0; i_site < sitenum; ++i_site) {
        for (int i_sat = 0; i_sat < sats[i_site].nsats_; ++ i_sat) {
            bindEph(sats[i_site].sat_[i_sat].sys_);
            Sattime tobs = sats[i_site].sat_[i_sat].obs_->time;
            double P = searchpseu(sats[i_site].sat_[i_sat].obs_->P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site].sat_[i_sat].sys_ << "prn: " << 
                        sats[i_site].sat_[i_sat].prn_ << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site].sat_[i_sat].clk[0];
            satpos_->satvel(tobs, sats[i_site].sat_[i_sat]);
            delete satpos_; satpos_ = nullptr;
        }
    }
    return true;
}

bool CPntbase::relativeeffect(sat* sats) {
    int sitenum = opt_->sitenum_;
    for (int i_site = 0; i_site < sitenum; ++ i_site)
        for (int i_sat = 0; i_sat < sats[i_site].nsats_; ++ i_sat) {
            bindEph(sats[i_site].sat_[i_sat].sys_);
            Sattime tobs = sats[i_site].sat_[i_sat].obs_->time;
            double P = searchpseu(sats[i_site].sat_[i_sat].obs_->P);
            if (P == 0)
                cout << "Warning, no pseudorange found, sys: " <<
                        sats[i_site].sat_[i_sat].sys_ << "prn: " << 
                        sats[i_site].sat_[i_sat].prn_ << endl;
            tobs = tobs - P / VEL_LIGHT - sats[i_site].sat_[i_sat].clk[0];
            satpos_->relfix(tobs, sats[i_site].sat_[i_sat]);
            delete satpos_; satpos_ = nullptr;
        }
}

void CPntbase::earthRotateFix(sat* sats) {
    int sitenum = opt_->sitenum_;
    for (int i_site = 0; i_site < sitenum; ++ i_site)
        for (int i_sat = 0; i_sat < sats[i_site].nsats_; ++ i_sat) {
            if (!sats[i_site].sat_[i_sat].eph_) continue;
            bindEph(sats[i_site].sat_[i_sat].sys_);
            Sattime dt = sats[i_site].sat_[i_sat].obs_->time -sats[i_site].sat_[i_sat].eph_->sig_;
            satpos_->earthRotateFix(dt._2sec(), sats[i_site].sat_[i_sat]);
            delete satpos_; satpos_  = nullptr;
        }
}

void CPntbase::satazel(double* site, sat &sat) {
    int satnum = sat.nsats_;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (*site == 0) {
            sat.sat_[i_sat].elev_ = sat.sat_[i_sat].azi_ = 0.0;
            break;
        } else {
            double neu[3];
            bindEph(sat.sat_[i_sat].sys_);
            XYZ2NEU(site, sat.sat_[i_sat].pos_, satpos_->GetElli(), neu);
            sat.sat_[i_sat].elev_ = atan(neu[2] / sqrt(neu[0] * neu[0] + neu[1] * neu[1]));
            sat.sat_[i_sat].azi_ = atan2(neu[0], neu[1]);
            delete satpos_; satpos_ = nullptr;
        }
    }
}

int CPntbase::usesys(sat &sats) {
    int satnum = sats.nsats_, num = 0;
    int used = 0x00;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat)
        used |= sats.sat_[i_sat].sys_;
    for (int i = 0; i < MAXSYS; ++i) 
        if ((used & SYS_ARRAY[i]) == SYS_ARRAY[i])
            num ++;
    return num;
}

void CPntbase::GetLC(sat sat, MatrixXd &w) {
    int satnum = sat.nsats_, num = 0;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sat.sat_[i_sat].isused) continue;
        if (!sat.sat_[i_sat].isused) continue;
        int p_pos[2] = {0};  // used pseudorange position
        for (int i = 0, j = 0; i < MAXFREQ; ++ i) 
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i])
                p_pos[j++] = i;
        double L1 = GetFreq(sat.sat_[i_sat].sys_, FREQ_ARRAY[p_pos[0]]);
        double L2 = GetFreq(sat.sat_[i_sat].sys_, FREQ_ARRAY[p_pos[1]]);
        double L1_square = L1 * L1, L2_square = L2 * L2;
        double IF = L1_square / (L1_square - L2_square) * sat.sat_[i_sat].obs_->P[p_pos[0]] -
                    L2_square / (L1_square - L2_square) * sat.sat_[i_sat].obs_->P[p_pos[1]];
        if(sat.sat_[i_sat].sys_ == SYS_BDS)
            IF -= VEL_LIGHT * L1_square * sat.sat_[i_sat].eph_->Tgd_[0] / 
                  (L1_square - L2_square);
        w(num , 0) += IF; num ++;
    }
}

void CPntbase::GetNonCombine(sat sat, MatrixXd &w) {
    int satnum = sat.nsats_, num = 0;
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sat.sat_[i_sat].isused) continue;
        int p_pos;  // used pseudorange position
        for (int i = 0; i < MAXFREQ; ++ i) 
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i])
                p_pos = i;
        double L1 = sat.sat_[i_sat].obs_->P[p_pos];
        w(num, 0) += L1; num ++;
    }
}


void CPntbase::setopt(prcopt opt) {
    // memcpy(&opt_, &opt, sizeof(opt));
    if (opt_) {
        opt_ = &opt;
    } else {
        memcpy(opt_, &opt, sizeof(opt));
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
    int satnum = sats.nsats_, num = 0;
    for (int isat = 0; isat < satnum; ++ isat) {
        if (!sats.sat_[isat].isused) continue;
        double T = Hopefield(sats.sat_[isat].elev_, H);
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
        if (res_->bpos_ecef_[0] != 0) {
            *xzy = res_->bpos_ecef_;
            *blh = res_->bpos_blh_;
        } else if (opt_->base_[0] != 0) {
            memcpy(res_->bpos_ecef_, opt_->base_, sizeof(double) * 3);
            *xzy = res_->bpos_ecef_;
            *blh = res_->bpos_blh_;
        } else {
            memset(res_->bpos_ecef_, 0, sizeof(double) * 3);
            *xzy = res_->bpos_ecef_;
            *blh = res_->bpos_blh_;
        }
    } else {
        if (res_->rpos_ecef_[0] != 0) {
            *xzy = res_->rpos_ecef_;
            *blh = res_->rpos_blh_;
        } else if (opt_->rover_[0] != 0) {
            memcpy(res_->rpos_ecef_, opt_->rover_, sizeof(double) * 3);
            *xzy = res_->rpos_ecef_;
            *blh = res_->rpos_blh_;
        } else {
            memset(res_->rpos_ecef_, 0, sizeof(double) * 3);
            *xzy = res_->rpos_ecef_;
            *blh = res_->rpos_blh_;
        }
    }
}

void CPntbase::getSysObs(sat sats, int* nobs) {
    if(!nobs) return;
    int nobs_sys[MAXSYS] = {0};
    for (int isat = 0; isat < sats.nsats_; ++isat) 
        for (int i = 0; i < MAXSYS; ++i) 
            if (sats.sat_[isat].sys_ == SYS_ARRAY[i] && sats.sat_[isat].isused)
                nobs_sys[i] += 1;
    memcpy(nobs, nobs_sys, sizeof(int) * MAXSYS);
}

double CPntbase::getSatUserPos(sat_s sats, double* sitepos) {
    double coef[3] = {0}, dist = 0;
    for (int i = 0; i < 3; ++i) {
        coef[i] = sitepos[i] - sats.pos_[i];
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
    res_->ambi_flag_ = FLOAT_SOLU;
    if (opt_->soltype_ == SOLTYPE_FLOAT) return;
    auto state = optimizer_->GetState();
    auto var = optimizer_->GetVar();
    // double sigma = optimizer_->GetInternalSigma(); // var = var * sigma;

    int n_ambi = state.row() - 3;
    auto float_solu = state.block<3, 0>(n_ambi, 1);
    auto float_solu_var = var.block<3, 3>(n_ambi, n_ambi);
    double* fix_tmp = new double[n_ambi * solu_num];
    MatrixXd fix(n_ambi, 1);
    MatrixXd s(1, solu_num);
    lambda(n_ambi, solu_num, float_solu._2array(), float_solu_var._2array(), fix_tmp, s._2array());
    for (int i = 0; i < n_ambi; ++ i) fix(i, 0) = fix_tmp[i];
    res_->ratio_ = s(0, 1) / s(0, 0);
    if (res_->ratio_ < RATIO_THRES)
        return;

    res_->ambi_flag_ = FIX_SOLU;
    int is_inverse = 0;
    // update result
    auto Qba = var.block<0, 3> (3, n_ambi);
    auto Qab = var.block<3, 0>(n_ambi, 3);
    auto Qbb = var.block<0, 0> (3, 3);
    auto state_float = state.block<0, 0>(3, 1);
    MatrixXd state_fix = state_float - Qba * float_solu_var.inverse(is_inverse) * (float_solu - fix);
    for (int i = 0; i < 3; ++ i) res_->rpos_ecef_[i] = state_fix(i, 0);
    delete[] fix_tmp; fix_tmp = nullptr;
}
