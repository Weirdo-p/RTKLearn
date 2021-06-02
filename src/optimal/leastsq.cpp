#include "navigation/optimal/leastsq.h"
#include "navigation/pnt/pntbase.h"
#include "navigation/pnt/pntrtk.h"

bool CLeastsq::optimize(MatrixXd B, MatrixXd P, MatrixXd w) {
    int istrue = 0;
    MatrixXd inv = B.transpose() * P * B;
    inv = inv.inverse(istrue);
    x_ = inv * B.transpose() * P * w;
    return bool(istrue);
}

MatrixXd CLeastsq::Getx() {
    return x_;
}

CLeastsq::CLeastsq() {

}

CLeastsq::CLeastsq(prcopt* opt) {
    if (!opt_) opt_ = opt;
}

void CLeastsq::getl(sat* sats, double* basepos, int* refsats, MatrixXd pos, MatrixXd &w) {
    int nobs = sats->nsats_;
    int ref_prn = -1, sysflag = -1;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
    #pragma omp parallel for
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
        double dist_base_ref = CPntbase::getSatUserPos(ref_sat_b, basepos);
        double dist_base_com = CPntbase::getSatUserPos(sats[0].sat_[isat], basepos);
        double dist_rover_ref = CPntbase::getSatUserPos(ref_sat_r, rover_pos);
        double dist_rover_com = CPntbase::getSatUserPos(sats[1].sat_[isat], rover_pos);
        double dist_cal = dist_rover_com - dist_rover_ref - dist_base_com + dist_base_ref;
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sats[1].sat_[isat].obs_->L[i] - ref_sat_r.obs_->L[i]-
                    sats[0].sat_[isat].obs_->L[i] + ref_sat_b.obs_->L[i];
                double pseu_obs = sats[1].sat_[isat].obs_->P[i] - ref_sat_r.obs_->P[i]- 
                    sats[0].sat_[isat].obs_->P[i] + ref_sat_b.obs_->P[i];
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
            for (int i = 0; i < 3; ++i) pos(i, 0) = res.rpos_ecef_[i];
        }
        double rover_pos[3] = {0};
        
        for (int i = 0; i < 3; ++i) rover_pos[i] = pos(i, 0);
        B.Zero(); P.Zero(); w.Zero();
        getDesign(sats_epoch, nobs, rover_pos, refsats, B);
        getl(sats_epoch, rover_pos, refsats, pos, w);
        getweight(sats_epoch, refsats, nobs, P);
        int flag = 0;
        if(!optimize(B, P, w)) {
            continue;
        }
        MatrixXd x = Getx();
        MatrixXd x_pos(3, 1);
        for (int i = 0; i < 3; ++i) {
            pos(i, 0) += x(i, 0); x_pos(i, 0) = x(i, 0);
        }
        if (x_pos.norm() < 1E-6) break;
        iter ++;
    }
    for (int i = 0; i < 3; ++ i) res.rpos_ecef_[i] = pos(i, 0);
}
