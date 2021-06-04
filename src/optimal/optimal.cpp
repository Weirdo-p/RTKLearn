#include "navigation/optimal/optimal.h"
#include "navigation/pnt/pntbase.h"

MatrixXd COptimal::getDoubleDiffCov(MatrixXd cov_sd, int* refpos, int* obs_sys, int nobs) {
    int nfreq = opt_->freqnum_;
    MatrixXd coeff((nobs - opt_->nsys_) * nfreq * 2, cov_sd.row());
    MatrixXd sub_coef(nfreq * 2, nfreq * 2);
    coeff.Zero(); sub_coef.Identity();
    int i_row = 0, col = 0;
    #pragma omp parallel for
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
        if (refpos[i_sys] >= total) col += 2 * nfreq;
    }

    return coeff * cov_sd * coeff.transpose();
}

MatrixXd COptimal::getSingleDiffCov(MatrixXd cov) {
    MatrixXd coeff(cov.row() / 2, cov.col());
    coeff.Zero();
    int nfreq = opt_->freqnum_;
    MatrixXd subcoef(2 * nfreq, 2 * nfreq);
    subcoef.Identity();
    int i_col = 0;
    // #pragma omp parallel for
    for (int i_row = 0; i_row < coeff.row(); i_row += subcoef.row()) {
        int count = 1;
        for (; i_col < coeff.col() && count <= 2; i_col += subcoef.col(), ++count) {
            coeff.block(i_row, i_col, subcoef * pow(-1, count));
        }
    }
    MatrixXd cov_diff = coeff * cov * coeff.transpose();
    return cov_diff;
}

void COptimal::getsubweight(sat_s sat_, MatrixXd &subcov) {
    subcov.resize(opt_->freqnum_ * 2, opt_->freqnum_ * 2);
    subcov.Zero();
    int num = 0;
    for (int i = 0; i < MAXFREQ; ++i) {
        if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]){
            subcov(num, num) = CPntbase::weightbyelev(sat_.elev_, 0.3, 1.5); num += 1;
            subcov(num, num) = CPntbase::weightbyelev(sat_.elev_, 0.001, 1.5); num += 1;
        }
    }
}

void COptimal::getweight(sat* sats, int* refsats, int nobs, MatrixXd &P) {
    int num = 0, ref_prn = 0, sysflag = 0;
    int nfreq = opt_->freqnum_, nsites = opt_->sitenum_;
    int refpos[MAXSYS] = {0}, obs_sys[MAXSYS] = {0};
    getSysObs(sats[1], obs_sys);
    MatrixXd cov(nobs * 2 * nfreq * nsites, nobs * 2 * nfreq * nsites);
    cov.Zero();
    #pragma omp parallel for
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

    int istrue = true;
    P = cov_double_diff.inverse(istrue);
}

int COptimal::findFreqPos(int sysflag, int* obs_sys, int &last_freq, double &freq) {
    int extern_pos = 0;
    if (opt_->freqnum_ == 1)
        extern_pos += 0;
    else {
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i] && 
                    (last_freq == -1 || last_freq != FREQ_ARRAY[i])) {
                last_freq = FREQ_ARRAY[i];
                freq = CPntbase::GetFreq(sysflag, last_freq);
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

int COptimal::findSysPos(int sysflag, int* obs_sys) {
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

void COptimal::getSysObs(sat sats, int* nobs) {
    if(!nobs) return;
    int nobs_sys[MAXSYS] = {0};
    for (int isat = 0; isat < sats.nsats_; ++isat) 
        for (int i = 0; i < MAXSYS; ++i) 
            if (sats.sat_[isat].sys_ == SYS_ARRAY[i] && sats.sat_[isat].isused)
                nobs_sys[i] += 1;
    memcpy(nobs, nobs_sys, sizeof(int) * MAXSYS);
}

void COptimal::getDesign(sat* sats, int nobs, double* sitepos, int* refsats, MatrixXd &B) {
    int num = 0, ref_prn = 0, sysflag = 0, pos = 0;
    int obs_sys[MAXSYS] = {0};
    getSysObs(sats[1], obs_sys);
    #pragma omp parallel for
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
        // #pragma omp parallel for
        for (int i_freq = 0; i_freq < MAXFREQ; ++ i_freq) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i_freq]) != FREQ_ARRAY[i_freq])
                continue;
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
}

sat_s COptimal::findRef(sat sats, int sysflag, int prn) {
    int nobs = sats.nsats_;
    for (int isat = 0; isat < nobs; ++isat) {
        if (sats.sat_[isat].sys_ == sysflag && sats.sat_[isat].prn_ == prn)
            return sats.sat_[isat];
    }
}

void COptimal::getDesignDim(sat sats, int nobs, int &row, int &col) {
    int nsys = opt_->nsys_, nfreq = opt_->freqnum_;
    row = col = 0;    
    row = (nobs - opt_->nsys_) * nfreq * 2;
    col = (nobs - opt_->nsys_) * nfreq + 3;
}

int COptimal::chooseref(sat sats, int* refsat) {
    static int REFSATS[MAXSYS] = { -1 };
    static int count = 0;
    count += 1;
    for (int i = 0; i < MAXSYS; ++i){
        if ((opt_->navsys_ & SYS_ARRAY[i]) == SYS_ARRAY[i])
            refsat[i] = chooseref(sats, SYS_ARRAY[i]);
    }
#if 0
    refsat[0] = 2; refsat[1] = 1;// for test
    if (count == 500)
        refsat[1] = 22;
#endif
    memcpy(REFSATS, refsat, sizeof(int) * MAXSYS);
    return 1;
}

int COptimal::chooseref(sat sats, int sysflag) {
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

COptimal::COptimal() {
    opt_ = nullptr;
}

void COptimal::setopt(prcopt* opt) {
    if (!opt) return;

    if (!opt_) {
        opt_ = opt;
    }
}

COptimal::COptimal(prcopt* opt) {
    if (!opt_) 
        opt_ = opt;
}

COptimal::~COptimal() {
    if (opt_)
        delete opt_;
    opt_ = nullptr;
}