#include "navigation/pnt/pntspp.h"
#include "navigation/atmosphere.h"
#include "navigation/coors.h"
#include "navigation/timemodule.h"
#include <iomanip>
#include <fstream>

CPntspp::CPntspp() {

}

CPntspp::CPntspp(prcopt opt) {
    if (opt_)
        memcpy(opt_, &opt, sizeof(opt));
    else {
        opt_ = new prcopt;
        memcpy(opt_, &opt, sizeof(opt));
    }
}

CPntspp::~CPntspp() {
    
}


int CPntspp::process() {
    
}

int CPntspp::spp(sat* sats) {
    int sitenum = opt_->sitenum_;
    int start = 0;
    if (opt_->base_[0] != 0) start = 1;
    for (int i_site = start; i_site < opt_->sitenum_; ++ i_site) {
        spp_site(i_site, &sats[i_site]);
    }    
}

int CPntspp::spp_site(int i_site, sat* sats) {
    int MAXITER = 10;
    double* sitepos_ecef, *sitepos_blh;

    selectPos(&sitepos_ecef, &sitepos_blh, i_site);
    int nobs = excludesats(*sats);
    double H;
    MatrixXd pos;
    Ellipsoid type(WGS84);
    for (int iter = 0; iter < MAXITER; ++ iter) {
        satazel(sitepos_ecef, *sats);
        if (i_site == 0) {
            H = res_->bpos_blh_[2];
        } else {
            H = res_->rpos_blh_[2];
        }
        int nsys = usesys(*sats);
        if (nsys != opt_->nsys_)
            cout << "WARNING: only " << nsys << " system(s) are used" << endl;
        int cols = 4;
        for (int i = 0; i < nsys - 1; ++i)
            cols += 1;
        if (nobs < cols) {
            cout << "FATAL: no sufficient observations" << endl;
            continue;
        }
        if (iter == 0) {
            pos.resize(cols, 1); pos.Zero();
            for (int i = 0; i < 3; ++i) pos(i, 0) = sitepos_ecef[i];
        }
        MatrixXd B(nobs, cols), P(nobs, nobs);
        P.Identity(); B.Zero();
        MatrixXd w(nobs, 1);
        Getl(i_site, *sats, sitepos_ecef, pos, w);
        // cout << w << endl << endl;
        GetDesign(*sats, sitepos_ecef, B);
        // cout << B << endl << endl;
        // cout << pos << endl << endl;
        CLeastsq optimizer(opt_);
        if(!optimizer.optimize(B, P, w, nobs)) {
            memset(sitepos_ecef, 0, sizeof(double) * 3);
            cout << w << endl;
            cout << B << endl;
            continue;
        }
        MatrixXd x = optimizer.Getx();
        for (int i = 0; i < x.row(); ++i) pos(i, 0) += x(i, 0);
        setres(sitepos_ecef, sitepos_blh, type, pos);
        if (x.norm() < 1E-6)
            break;
    }
    sitepos_ecef = nullptr; sitepos_blh = nullptr;
}

void CPntspp::GetDesign(sat sats, double* sitepos, MatrixXd &B) {
    int satnum = sats.nsats_, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats.sat_[i_sat].isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats.sat_[i_sat].pos_[i];
            dist += coeff[i] * coeff[i];
        }
        int p_pos = 0;  // used pseudorange position
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats.sat_[i_sat].sys_ & SYS_ARRAY[i]) == SYS_ARRAY[i])
                p_pos = i;
        if (opt_->nsys_ == 1) p_pos = 0;
        dist = sqrt(dist);
        for (int i = 0; i < B.col(); ++i) {
            if (i < 3)
                B(num, i) = coeff[i] / dist;
            else {
                B(num, i + p_pos) = 1;
                break;
            }
        }
        num++;
    }
}

void CPntspp::Getl(int i_site, sat sats, double* sitepos, MatrixXd pos, MatrixXd &w) {
    int satnum = sats.nsats_, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats.sat_[i_sat].isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats.sat_[i_sat].pos_[i];
            dist += coeff[i] * coeff[i];
        }
        int p_pos = 0;
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats.sat_[i_sat].sys_ & SYS_ARRAY[i]) == SYS_ARRAY[i])
                p_pos = i;
        if (opt_->nsys_ == 1) p_pos = 0;
        dist = sqrt(dist);
        w(num++, 0) = -dist + sats.sat_[i_sat].clk[0] * VEL_LIGHT - pos(3 + p_pos, 0);
    }
    if (opt_->freqnum_ != 1) GetLC(sats, w);
    else GetNonCombine(sats, w);
    double BLH[3];
    if (*sitepos == 0) return;
    XYZ2BLH(sitepos, Ellipsoid(WGS84), BLH);
    if (abs(BLH[2]) > 1E4) 
        return;
    tropfix(sats, w, BLH[2]);
}

res_t* CPntspp::getRes() {
    return res_;
}

int CPntspp::excludesats(sat &sat) {
    int satnum = sat.nsats_;
    int n = 0;
    double cutoff = opt_->elecutoff_;
    for(int isat = 0; isat < satnum; ++ isat) {
        if(sat.sat_[isat].elev_ != 0 && sat.sat_[isat].elev_ < opt_->elecutoff_) {
            sat.sat_[isat].isused = false; continue;
        }
        if (!sat.sat_[isat].eph_) {
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
                n += 1;
                sat.sat_[isat].isused = true;
            }
            else sat.sat_[isat].isused = false;
        }
    }

    return n;
}
