#include "navigation/pnt/pntrtk.h"
#include "navigation/timemodule.h"
#include "navigation/coors.h"
#include <string.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <chrono>
#include "navigation/optimal/rtkekf.h"
using namespace chrono;

CPntrtk::CPntrtk() {
    opt_ = nullptr; optimizer_ = nullptr;
}

CPntrtk::CPntrtk(prcopt opt) {
    if (opt_) {
        memcpy(opt_, &opt, sizeof(prcopt));
    } else {
        opt_ = new prcopt;
        memcpy(opt_, &opt, sizeof(prcopt));
    }
    if (!optimizer_) {
        if (opt_->proctype_ == PROC_KF)
            optimizer_ = new CRtkekf(opt_);
        if (opt_->proctype_ == PROC_LS)
            optimizer_ = new CLeastsq(opt_);
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
        spprunner_->spp(sats_epoch);
        auto res = spprunner_->getRes();
        memcpy(res_->rpos_ecef_, res->rpos_ecef_, sizeof(double) * 3);
        satazel(res_->rpos_ecef_, sats_epoch[0]);
        rtk(sats_epoch);
        auto end = system_clock::now();
        auto cost = duration_cast<microseconds> (end - start);
        cout << setw(4) << sats_epoch->sat_->obs_->time.Week_ << " " << setw(9) << fixed << setprecision(6) << sats_epoch->sat_->obs_->time.Sow_;
        cout<< " TIME COST: " << double (cost.count()) * microseconds::period::num / microseconds::period::den << " seconds" << endl;
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
    #pragma omp parallel for
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
        #pragma omp parallel for
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
    #pragma omp parallel for
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

bool CPntrtk::rtk(sat* sats_epoch) {
    int MAXITER = 10, iter = 0, issuccess = 0;
    int REFSATS[MAXSYS] = { -1 };
    int nobs = 0;
    for (int i = 0; i < opt_->sitenum_; ++i)
        excludesats(sats_epoch[i]);  // observation count by rover
    
    nobs = obsnumber(sats_epoch);

    optimizer_->optimize(sats_epoch, *res_);
    XYZ2NEU(res_->bpos_ecef_, res_->rpos_ecef_, WGS84, res_->enu);
    ofstream out("./out1.txt", ios::app);
    out << setfill(' ');
    out << setw(4) << sats_epoch->sat_->obs_->time.Week_ << " " << setw(18) << fixed << setprecision(6) << sats_epoch->sat_->obs_->time.Sow_;
    out << " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[0] <<
            " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[1] <<
            " " << setw(18) << fixed << setprecision(6) << res_->rpos_ecef_[2] << 
            " " << setw(9) << fixed << setprecision(6) << res_->enu[0] <<
            " " << setw(9) << fixed << setprecision(6) << res_->enu[1] <<
            " " << setw(9) << fixed << setprecision(6) << res_->enu[2] << " ";
    out << setw(3) << fixed << nobs << endl;
    out.close();
}

int CPntrtk::obsnumber(sat* sats) {
    int nobs = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < sats->nsats_; ++ i_sat) {
        if (!sats[0].sat_[i_sat].isused || !sats[1].sat_[i_sat].isused)
            sats[0].sat_[i_sat].isused = sats[1].sat_[i_sat].isused = false;
        else nobs += 1;
    }
    return nobs;
}

int CPntrtk::excludesats(sat &sat) {
    int satnum = sat.nsats_;
    double cutoff = opt_->elecutoff_;
    for(int isat = 0; isat < satnum; ++ isat) {
        if(sat.sat_[isat].elev_ != 0 && sat.sat_[isat].elev_ < opt_->elecutoff_) {
            sat.sat_[isat].isused = false; continue;
        }
        if ((sat.sat_[isat].prn_ <= 5) && sat.sat_[isat].sys_ == SYS_BDS) {
            sat.sat_[isat].isused = false; continue;
        }
        if ((sat.sat_[isat].prn_ >= 59) && sat.sat_[isat].sys_ == SYS_BDS) {
            sat.sat_[isat].isused = false; continue;
        }

        if(sat.sat_[isat].prn_ == 5 && sat.sat_[isat].sys_ == SYS_BDS) {
            sat.sat_[isat].isused = false; continue;
        }

        else {
            bool isobs = true;
            for (int i = 0; i < MAXFREQ; ++i)
                if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i] && 
                    abs(sat.sat_[isat].obs_->P[i]) <= 1e-6) {
                        isobs = false; break;
                    }
            if (isobs) {
                sat.sat_[isat].isused = true;
            }
            else sat.sat_[isat].isused = false;
        }

    }
}
