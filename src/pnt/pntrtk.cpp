#include "navigation/pnt/pntrtk.h"
#include "navigation/timemodule.h"
#include <string.h>
#include <stdio.h>
#include <iomanip>

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
    double REFSATS[MAXSYS] = { -1 };
    if (opt_->sitenum_ < 2) return -1;
    sat* sats_epoch = new sat[opt_->sitenum_];
    while(inputobs(sats_epoch)) {
        satclk(sats_epoch);
        satpos(sats_epoch);
        earthRotateFix(sats_epoch);
        satvel(sats_epoch);
        relativeeffect(sats_epoch);
        if (res_->bpos_ecef_[0] == 0 || res_->rpos_ecef_[0] == 0) {
            spprunner_->spp(sats_epoch);
            res_ = spprunner_->getRes();
        }
        chooseref(sats_epoch[1], REFSATS);
        int nobs = excludesats(sats_epoch[1]);  // observation count by rover
        int row = 0, col = 0;
        getDesignDim(nobs, row, col);
        MatrixXd B(row, col);
        MatrixXd w(row, 1);
        MatrixXd P(row, row);
        B.Zero(); P.Zero(); w.Zero();
        getDesign(sats_epoch, nobs, res_->rpos_ecef_, REFSATS, B);
        memset(sats_epoch, 0, sizeof(sat) * opt_->sitenum_);
    }
    cout << "test" << endl;
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

void CPntrtk::getDesignDim(int nobs, int &row, int &col) {
    int nsys = opt_->nsys_, nfreq = opt_->freqnum_;
    row = col = 0;
    int basic_row = nobs - 1;    // basic dimension
    row = basic_row * nsys * nfreq;
    col = basic_row * nsys * nfreq + 3;
}

void CPntrtk::getDesign(sat* sats, int nobs, double* sitepos, double* refsats, MatrixXd &B) {
    int num = 0, ref_prn = 0, sysflag = 0;
    for (int isat = 0; isat < sats->nsats_; ++ isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((opt_->navsys_ & SYS_ARRAY[i]) == SYS_ARRAY[i]) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1].sat_[isat].prn_ == ref_prn) continue;
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
        // fill in B with pseudorange and phase observations
        for(int count = 0; count < 2; ++count, ++num)
            for (int i = 0; i < 3; ++i)
                B(num, i) = -coeff_sat[i] / dist_sat + coeff_ref[i] / dist_ref;
        // left part of B
    }
}

sat_s CPntrtk::findRef(sat sats, int sysflag, int prn) {
    int nobs = sats.nsats_;
    for (int isat = 0; isat < nobs; ++isat) {
        if (sats.sat_[isat].sys_ == sysflag && sats.sat_[isat].prn_ == prn)
            return sats.sat_[isat];
    }
}
