#include "navigation/pnt/pntrtk.h"
#include "navigation/timemodule.h"
#include "navigation/coors.h"
#include <string.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <chrono>
using namespace chrono;

CPntrtk::CPntrtk() {

}

CPntrtk::CPntrtk(prcopt opt) {
    if (opt_)
        memcpy(opt_, &opt, sizeof(opt));
    else {
        opt_ = new prcopt;
        memcpy(opt_, &opt, sizeof(opt));
    }
    spprunner_ = new CPntspp(*opt_);
}

CPntrtk::~CPntrtk() {
    if (spprunner_) 
        delete spprunner_;
    spprunner_ = nullptr;
}


int CPntrtk::process() {
    if (opt_->sitenum_ < 2) return -1;
    sat* sats_epoch = new sat[opt_->sitenum_];
    if (opt_->base_[0] != 0) memcpy(res_->bpos_ecef_, opt_->base_, sizeof(double) * 3);
    if (opt_->rover_[0] != 0) memcpy(res_->rpos_ecef_, opt_->rover_, sizeof(double) * 3);
    while(inputobs(sats_epoch)) {
        auto start = system_clock::now();
        satclk(sats_epoch);
        satpos(sats_epoch);
        earthRotateFix(sats_epoch);
        satvel(sats_epoch);
        relativeeffect(sats_epoch);
        spprunner_->spp(&sats_epoch);
        auto res = spprunner_->getRes();
        memcpy(res_->rpos_ecef_, res->rpos_ecef_, sizeof(double) * 3);
        satazel(res_->rpos_ecef_, sats_epoch[0]);
        rtk(sats_epoch);
        auto end = system_clock::now();
        auto cost = duration_cast<microseconds> (end - start);
        cout << setw(4) << sats_epoch->sat_->obs_->time.Week_ << " " << setw(18) << fixed << setprecision(6) << sats_epoch->sat_->obs_->time.Sow_;
        cout<< " TIME COST: " << double (cost.count()) * microseconds::period::num / microseconds::period::den << "seconds" << endl;
        memset(sats_epoch, 0, sizeof(sat) * opt_->sitenum_);
    }
}

int CPntrtk::inputobs(sat* sats) {
    // record observations' postion
    static int base_pos = 0, rover_pos = 0;
    int nsites = opt_->sitenum_;
    obs* obss = rnx_->GetObs();
    nav* navs = rnx_->GetEph();
    obs* base = &obss[0], *rover = &obss[1];
    double diff = Sattimediff(base->obs_[base_pos].time, rover->obs_[rover_pos].time);
    // searching observations
    while(abs(diff) >= 1e-3) {
        if (base_pos >= base->obsnum_ || rover_pos >= rover->obsnum_)
            break;
        if (diff < 0)   base_pos += 1;
        else            rover_pos += 1;
        diff = Sattimediff(base->obs_[base_pos].time, rover->obs_[rover_pos].time);
    }
    // searching ephemeris
    if (base_pos >= base->obsnum_ || rover_pos >= rover->obsnum_)
        return 0;
    setobs(obss, base_pos, rover_pos, sats);
    setnav(navs, sats);
    obss = nullptr; navs = nullptr; base = nullptr; rover = nullptr;
    return 1;
}

int CPntrtk::setobs(obs* obss, int &base_pos, int &rover_pos, sat* sats) {
    obs* base = &obss[0], *rover = &obss[1];
    int nsats = 0;
    Sattime time = base->obs_[base_pos].time;
    int i_base = base_pos;
    for (; i_base < base->obsnum_; ++ i_base) {
        if (base->obs_[i_base].time != time) break;
        for (int i_rover = rover_pos; i_rover < rover->obsnum_; ++ i_rover) {
            if (rover->obs_[i_rover].time != time) break;
            if (base->obs_[i_base].sys != rover->obs_[i_rover].sys || 
                base->obs_[i_base].sat != rover->obs_[i_rover].sat)
                continue;
            sats[0].sat_[nsats].prn_ = sats[1].sat_[nsats].prn_ = rover->obs_[i_rover].sat;
            sats[0].sat_[nsats].sys_ = sats[1].sat_[nsats].sys_ = rover->obs_[i_rover].sys;
            sats[0].sat_[nsats].obs_ = &base->obs_[i_base];
            sats[1].sat_[nsats].obs_ = &rover->obs_[i_rover];
            nsats += 1;
        }
    }
    sats[0].nsats_ = sats[1].nsats_ = nsats;
    base_pos = i_base;
    base = nullptr; rover = nullptr;
    return 1;
}

int CPntrtk::setnav(nav* navs, sat* sats) {
    for (int i = 0; i < 2; ++i) {
        for (int isat = 0; isat < sats[i].nsats_; ++isat) {
            int sys = sats[i].sat_[isat].sys_;
            int prn = sats[i].sat_[isat].prn_;
            Sattime time = sats[i].sat_[isat].obs_->time;
            sats[i].sat_[isat].eph_ = searchnav(time, sys, prn, navs);
        }
    }
    return 1;
}

nav_t* CPntrtk::searchnav(Sattime time, int sys, int prn, nav* navs) {
    int delay = 0;
    bool isfind = false;
    if (sys == SYS_GPS) delay = 7200;
    else if (sys == SYS_BDS) {
        time = time - BDT2GPST;
        delay = 3600;
    }
    // nav_t* sat_nav = nullptr;
    for (int i_nav = navs->num; i_nav > 0; -- i_nav) {
        if (sys != navs->msg_[i_nav].sys_ || prn != navs->msg_[i_nav].prn_)
            continue;
        if (abs(Sattimediff(time, navs->msg_[i_nav].toe_)) > delay)
            continue;
        isfind = true;
        return &navs->msg_[i_nav];
    }
    if (!isfind) {
        cout << "no ephemeris found: sys: " << sys << " prn: " << prn << endl;
        return nullptr;
    }
}

int CPntrtk::chooseref(sat sats, int sysflag) {
    int satsnum = sats.nsats_, prn = 0;
    double max_elev = -1;
    for (int i_sat = 0; i_sat < satsnum; ++ i_sat) {
        if (sats.sat_[i_sat].isused &&
            sats.sat_[i_sat].elev_ >= max_elev && 
            sats.sat_[i_sat].sys_ == sysflag) {
            max_elev = sats.sat_[i_sat].elev_; prn = sats.sat_[i_sat].prn_;
        }
    }
    return prn;
}


int CPntrtk::chooseref(sat sats, double* refsat) {
    static double REFSATS[MAXSYS] = { -1 };
    for (int i = 0; i < MAXSYS; ++i){
        if ((opt_->navsys_ & SYS_ARRAY[i]) == SYS_ARRAY[i])
            refsat[i] = chooseref(sats, SYS_ARRAY[i]);
    }
    for (int i = 0; i < MAXSYS; ++i) 
        if (REFSATS[i] != refsat[i]) {
            // may be reset the kalman filter
        }
    memcpy(REFSATS, refsat, sizeof(double) * MAXSYS);
    return 1;
}

void CPntrtk::getDesignDim(sat sats, int nobs, int &row, int &col) {
    int nsys = opt_->nsys_, nfreq = opt_->freqnum_;
    row = col = 0;    
    row = (nobs - opt_->nsys_) * nfreq * 2;
    col = (nobs - opt_->nsys_) * nfreq + 3;
}

void CPntrtk::getDesign(sat* sats, int nobs, double* sitepos, double* refsats, MatrixXd &B) {
    int num = 0, ref_prn = 0, sysflag = 0, pos = 0;
    int obs_sys[MAXSYS] = {0};
    getSysObs(sats[1], obs_sys);
    for (int isat = 0; isat < sats->nsats_; ++ isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((opt_->navsys_ & SYS_ARRAY[i]) == sats[1].sat_[isat].sys_) {
                if (sysflag != SYS_ARRAY[i]) pos = 0;
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1].sat_[isat].prn_ == ref_prn && sats[1].sat_[isat].sys_ == sysflag) continue;
        sat_s ref_sat = findRef(sats[1], sysflag, ref_prn);
        double dist_sat = 0, dist_ref = 0, coeff_sat[3] = {0}, coeff_ref[3] = {0};
        for (int i = 0; i < 3; i ++) {
            // for reference satellite
            coeff_ref[i] = ref_sat.pos_[i] - sitepos[i];
            dist_ref += coeff_ref[i] * coeff_ref[i];
            // for the other one
            coeff_sat[i] = sats[1].sat_[isat].pos_[i] - sitepos[i];
            dist_sat += coeff_sat[i] * coeff_sat[i];
        }
        dist_sat = sqrt(dist_sat); dist_ref = sqrt(dist_ref);

        int last_freq = -1;
        for (int i_freq = 0; i_freq < opt_->freqnum_; ++ i_freq) {
            // fill in B with pseudorange and phase observations
            for(int count = 0; count < 2; ++count, ++num)
                for (int i = 0; i < 3; ++i)
                    B(num, i) = -coeff_sat[i] / dist_sat + coeff_ref[i] / dist_ref;
            // left part of B
            double freq = 0;
            int base_pos = 3, extern_pos = 0;
            extern_pos += findSysPos(sysflag, obs_sys);
            extern_pos += findFreqPos(sysflag, obs_sys, last_freq, freq);
            B(num - 1, base_pos + extern_pos + pos) = VEL_LIGHT / freq;
        }
        pos += 1;
    }
    // ofstream out("./B_debug.txt");
    // out << B << endl;
    // out.close();
}

sat_s CPntrtk::findRef(sat sats, int sysflag, int prn) {
    int nobs = sats.nsats_;
    for (int isat = 0; isat < nobs; ++isat) {
        if (sats.sat_[isat].sys_ == sysflag && sats.sat_[isat].prn_ == prn)
            return sats.sat_[isat];
    }
}

int CPntrtk::findSysPos(int sysflag, int* obs_sys) {
    int extern_pos = 0;
    if (opt_->nsys_ == 1)
        extern_pos += 0;
    else
        for (int i = 0; i < MAXSYS; ++ i) 
            if (SYS_ARRAY[i] == sysflag)
                if(i == 0)  extern_pos += 0;
                else
                    for (int j = i - 1; j >= 0; --j)
                        if (obs_sys[j] != 0)
                            extern_pos += (obs_sys[j] - 1) * opt_->freqnum_;
    return extern_pos;
}

int CPntrtk::findFreqPos(int sysflag, int* obs_sys, int &last_freq, double &freq) {
    int extern_pos = 0;
    if (opt_->freqnum_ == 1)
        extern_pos += 0;
    else {
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i] && 
                    (last_freq == -1 || last_freq != FREQ_ARRAY[i])) {
                last_freq = FREQ_ARRAY[i];
                freq = GetFreq(sysflag, last_freq);
                for (int isys = 0; isys < MAXSYS; ++ isys) 
                    if (SYS_ARRAY[isys] == sysflag)
                        if(i == 0)  extern_pos += 0;
                        else extern_pos += (obs_sys[isys] - 1);
                break;
            }
        }
    }
    return extern_pos;
}

void CPntrtk::getl(sat* sats, double* sitepos, double* refsats, MatrixXd pos, MatrixXd &w) {
    int nobs = sats->nsats_;
    int ref_prn = -1, sysflag = -1;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((opt_->navsys_ & SYS_ARRAY[i]) == sats[1].sat_[isat].sys_) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1].sat_[isat].prn_ == ref_prn && sats[1].sat_[isat].sys_ == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        double dist_base_ref = getSatUserPos(ref_sat_b, res_->bpos_ecef_);
        double dist_base_com = getSatUserPos(sats[0].sat_[isat], res_->bpos_ecef_);
        double dist_rover_ref = getSatUserPos(ref_sat_r, rover_pos);
        double dist_rover_com = getSatUserPos(sats[1].sat_[isat], rover_pos);
        double dist_cal = dist_rover_com - dist_rover_ref - dist_base_com + dist_base_ref;
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sats[1].sat_[isat].obs_->L[i] - ref_sat_r.obs_->L[i]-
                    sats[0].sat_[isat].obs_->L[i] + ref_sat_b.obs_->L[i];
                double pseu_obs = sats[1].sat_[isat].obs_->P[i] - ref_sat_r.obs_->P[i]- 
                    sats[0].sat_[isat].obs_->P[i] + ref_sat_b.obs_->P[i];
                double freq = GetFreq(sysflag, FREQ_ARRAY[i]);
                w(num ++, 0) = (pseu_obs - dist_cal);
                w(num ++, 0) = (phase_obs * (VEL_LIGHT / freq) - dist_cal);
            }
        }
    }
}

void CPntrtk::getweight(sat* sats, double* refsats, int nobs, MatrixXd &P) {
    int satnum = sats->nsats_;
    int num = 0, ref_prn = 0, sysflag = 0;
    int nfreq = opt_->freqnum_, nsites = opt_->sitenum_;
    int refpos[MAXSYS] = {0}, obs_sys[MAXSYS] = {0};
    getSysObs(sats[1], obs_sys);
    MatrixXd cov(nobs * 2 * nfreq * nsites, nobs * 2 * nfreq * nsites);
    cov.Zero();
    for (int isat = 0; isat < sats->nsats_; ++ isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) 
            if (refsats[i] == sats->sat_[isat].prn_ && SYS_ARRAY[i] == sats->sat_[isat].sys_)
                refpos[i] = num / 2;
        for (int isite = 0; isite < nsites; ++ isite) {
            MatrixXd subcov;
            getsubweight(sats[isite].sat_[isat], subcov);
            for (int i = 0; i < subcov.row(); ++ i, ++ num) 
                for (int j = 0; j < subcov.col(); ++ j)
                    if (i == j)
                        cov(num, num) = subcov(i, j);
        }
    }
    MatrixXd cov_single_diff = getSingleDiffCov(cov);
    MatrixXd cov_double_diff = getDoubleDiffCov(cov_single_diff, refpos, obs_sys, nobs);
    // cout << cov.row() << " " << cov.col() << " " << nobs << " ";
    // cout << cov_single_diff.row() << " " << cov_single_diff.col() << " " << cov_double_diff.row() << " " << cov_double_diff.col() << endl;

    int istrue = true;
    P = cov_double_diff.inverse(istrue);
    if (!istrue) {
        ofstream out ("./p_debug.txt", ios::app);
        out << P << endl << endl;
        out << cov_single_diff << endl << endl;
        out << cov_double_diff << endl << endl;
        out << cov << endl << endl;
        out.close();
    }
}

void CPntrtk::getsubweight(sat_s sat_, MatrixXd &subcov) {
    subcov.resize(opt_->freqnum_ * 2, opt_->freqnum_ * 2);
    subcov.Zero();
    int num = 0;
    for (int i = 0; i < MAXFREQ; ++i) {
        if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]){
            subcov(num, num) = weightbyelev(sat_.elev_, 0.01, 1.5); num += 1;
            subcov(num, num) = weightbyelev(sat_.elev_, 0.003, 1.5); num += 1;
        }
    }
}

MatrixXd CPntrtk::getSingleDiffCov(MatrixXd cov) {
    MatrixXd coeff(cov.row() / 2, cov.col());
    coeff.Zero();
    int nfreq = opt_->freqnum_;
    MatrixXd subcoef(2 * nfreq, 2 * nfreq);
    subcoef.Identity();
    int i_col = 0;
    for (int i_row = 0; i_row < coeff.row(); i_row += subcoef.row()) {
        int count = 1;
        for (; i_col < coeff.col() && count <= 2; i_col += subcoef.col(), ++count) {
            coeff.block(i_row, i_col, subcoef * pow(-1, count));
        }
    }
    MatrixXd cov_diff = coeff * cov * coeff.transpose();
    return cov_diff;
}

MatrixXd CPntrtk::getDoubleDiffCov(MatrixXd cov_sd, int* refpos, int* obs_sys, int nobs) {
    int nfreq = opt_->freqnum_;
    MatrixXd coeff((nobs - opt_->nsys_) * nfreq * 2, cov_sd.row());
    MatrixXd sub_coef(nfreq * 2, nfreq * 2);
    coeff.Zero(); sub_coef.Identity();
    int i_row = 0, col = 0;
    
    for (int i_sys = 0; i_sys < MAXSYS; ++i_sys) {
        int total = 0;
        for (int i = i_sys; i >= 0; i--) total += (obs_sys[i] - 1) * 2 * nfreq;
        for (; i_row < coeff.row() && i_row < total; i_row += 2 * nfreq, col += 2 * nfreq) {
            if (refpos[i_sys] == col) {
                i_row -= 2 * nfreq;
                continue;
            }
            coeff.block(i_row, refpos[i_sys], sub_coef * -1);
            coeff.block(i_row, col, sub_coef);
        }
    }
    // cout << nobs << " coeff" << " " << coeff.row() << " " << coeff.col() << "  " << (nobs - opt_->nsys_) * nfreq * 2 << endl;
    return coeff * cov_sd * coeff.transpose();
}

bool CPntrtk::rtk(sat* sats_epoch) {
    int MAXITER = 10, iter = 0;
    double REFSATS[MAXSYS] = { -1 };
    int nobs = 0;
    for (int i = 0; i < opt_->sitenum_; ++i)
        excludesats(sats_epoch[i]);  // observation count by rover
    nobs = obsnumber(sats_epoch);
    chooseref(sats_epoch[1], REFSATS);
    int row = 0, col = 0;
    MatrixXd pos;
    getDesignDim(sats_epoch[1], nobs, row, col);
    while (iter < MAXITER) {
        MatrixXd B(row, col);
        MatrixXd w(row, 1);
        MatrixXd P(row, row);
        if (iter == 0) {
            pos.resize(col, 1); pos.Zero();
            for (int i = 0; i < 3; ++i) pos(i, 0) = res_->rpos_ecef_[i];
        }
        double rover_pos[3] = {0};
        for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
        B.Zero(); P.Zero(); w.Zero();
        getDesign(sats_epoch, nobs, rover_pos, REFSATS, B);
        getl(sats_epoch, rover_pos, REFSATS, pos, w);
        getweight(sats_epoch, REFSATS, nobs, P);
        if(!optimizer_.optimize(B, P, w)) {
            continue;
        }
        MatrixXd x = optimizer_.Getx();
        MatrixXd x_pos(3, 1);
        for (int i = 0; i < 3; ++i) {
            pos(i, 0) += x(i, 0); x_pos(i, 0) = x(i, 0);
        }
        if (x_pos.norm() < 1E-8) break;
        iter ++;
    }
    for (int i = 0; i < 3; ++ i) res_->rpos_ecef_[i] = pos(i, 0);
    XYZ2NEU(res_->bpos_ecef_, res_->rpos_ecef_, WGS84, res_->enu);
    ofstream out("./out.txt", ios::app);
    out << setfill(' ');
    out << setw(4) << sats_epoch->sat_->obs_->time.Week_ << " " << setw(18) << fixed << setprecision(6) << sats_epoch->sat_->obs_->time.Sow_;
    out << " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[0] <<
            " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[1] <<
            " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[2] << 
            " " << setw(9) << fixed << setprecision(6) << res_->enu[0] <<
            " " << setw(9) << fixed << setprecision(6) << res_->enu[1] <<
            " " << setw(9) << fixed << setprecision(6) << res_->enu[2] << endl;
    out.close();
}

int CPntrtk::obsnumber(sat* sats) {
    int nsites = opt_->sitenum_;
    int nobs = 0;
    for (int i_sat = 0; i_sat < sats->nsats_; ++ i_sat) {
        if (!sats[0].sat_[i_sat].isused || !sats[1].sat_[i_sat].isused)
            sats[0].sat_[i_sat].isused = sats[1].sat_[i_sat].isused = false;
        else nobs += 1;
    }
    return nobs;
}
