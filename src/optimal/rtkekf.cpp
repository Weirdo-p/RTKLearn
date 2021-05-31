#include "navigation/optimal/rtkekf.h"

MatrixXd CRtkekf::filter(sat* sats_epoch, res_t &res, MatrixXd obs, int nobs, int* refsats, int* sysobs) {
    ofstream out("./kfdebug.txt");
    out_debug.open("./debug_before_fix.txt", ios::app); // TODO
    auto test = obssats(sats_epoch, nobs, refsats);
    out_debug << "EPOCH: " << sats_epoch->sat_->obs_->time.Week_ << "  " << sats_epoch->sat_->obs_->time.Sow_ << endl;
    out_debug << "observation numbers: " << nobs << endl;
    out_debug << "state dim: " << state_.row() << endl;
    out_debug << "last obssats: \n" << obssats_ << endl << endl;
    out_debug << "now obssats: \n" << test << endl << endl;
    out_debug.close();
    out << state_ << endl << endl;
    if (sats_epoch->sat_->obs_->time.Sow_ == 144090) {
        auto obs_sats_test = obssats(sats_epoch, nobs, refsats);
        out << obs_sats_test << endl << endl;
        cout << "test" << endl;
    }
    updateStatus(obs, sats_epoch, nobs, refsats, res, sysobs);
    statusFix(sats_epoch, nobs, refsats, sysobs);
    out_debug.open("./debug_after_fix.txt", ios::app); // TODO
    test = obssats(sats_epoch, nobs, refsats);
    out_debug << "EPOCH: " << sats_epoch->sat_->obs_->time.Week_ << "  " << sats_epoch->sat_->obs_->time.Sow_ << endl;
    out_debug << "observation numbers: " << nobs << endl;
    out_debug << "state dim: " << state_.row() << endl;
    out_debug << "last obssats: \n" << obssats_ << endl << endl;
    out_debug << "now obssats: \n" << test << endl << endl;
    out_debug.close();
    int issuccess = 0;
    MatrixXd state_time_predict = state_trans_ * state_;
    MatrixXd var_state_time_predict = state_trans_ * var_state_ * state_trans_.transpose() + var_sys_;
    MatrixXd approx(obs.row(), 1);
    getDDApprox(sats_epoch, state_time_predict, res.bpos_ecef_, refsats, approx);
    MatrixXd K = var_state_time_predict * design_.transpose() * (
        design_ * var_state_time_predict * design_.transpose() + var_obs_
    ).inverse(issuccess);
    state_ = state_time_predict + K * (obs - approx);
    MatrixXd temp = K * design_;
    MatrixXd identity(temp.row(), temp.col()); identity.Identity();
    var_state_ = (identity - temp) * var_state_time_predict;

    for (int i = 0; i < 3; ++ i)
        res.rpos_ecef_[i] = state_(i, 0);
    MatrixXd tmp = obs - approx;
    out << state_ << endl;
    out << tmp << endl << endl;
    out << test << endl << endl;
    out << obssats_ << endl << endl;
    // cout << setprecision(15) << state_ << endl;
    out.close();
}

void CRtkekf::getDDApprox(sat* sats, MatrixXd state, double* b_pos, int* refsats, MatrixXd &w) {
    int nobs = sats->nsats_;
    int ref_prn = -1, sysflag = -1;
    int freqnum = opt_->freqnum_;
    int num = 0, ar_pos_ = 3;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = state(i, 0);
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
        double dist_base_ref = dist(ref_sat_b.pos_, b_pos);
        double dist_base_com = dist(sats[0].sat_[isat].pos_, b_pos);
        double dist_rover_ref = dist(ref_sat_r.pos_, rover_pos);
        double dist_rover_com = dist(sats[1].sat_[isat].pos_, rover_pos);
        double dist_cal = dist_rover_com - dist_rover_ref - dist_base_com + dist_base_ref;
        int sys = 0;
        for (int i_freq = 0; i_freq < MAXFREQ; ++i_freq) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i_freq]) == FREQ_ARRAY[i_freq]) {
                double freq = GetFreq(sysflag, FREQ_ARRAY[i_freq]);
                w(num ++, 0) = dist_cal;
                // states are: pos, GPS_Nf1, GPS_Nf2, BDS_NF1, BDS_NF2
                for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                    int total = 0;
                    for (int i_t = i_sys; i_t >= 0; -- i_t)
                        total += (sysobs_[i_t] - 1) * 2 * freqnum;
                    if (SYS_ARRAY[i_sys] == sysflag) {
                        w(num ++, 0) = dist_cal + state(ar_pos_ + i_freq * (sysobs_[i_sys] - 1), 0) * (VEL_LIGHT / freq);
                        if (num == total)
                            ar_pos_ = total / 2 + 3 - 1;
                    }
                }
            }
        }
        ar_pos_ += 1;
    }
}

CRtkekf::CRtkekf() {
    memset(LeapSatsPos_, 0, sizeof(int) * MAXOBS);
    status_ = RTK_INIT;
    opt_ = nullptr;
}

MatrixXd CRtkekf::obssats(sat* sats_epoch, int nobs, int* refsats) {
    int nsats = sats_epoch->nsats_;
    MatrixXd usedsats(nobs, 3); usedsats.Zero();
    for (int i_sat = 0, i_row = 0; i_sat < nsats; ++i_sat) {
        if (!sats_epoch->sat_[i_sat].isused) continue;
        usedsats(i_row, 0) = sats_epoch->sat_[i_sat].sys_;
        usedsats(i_row, 1) = sats_epoch->sat_[i_sat].prn_;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) 
            if (SYS_ARRAY[i_sys] == sats_epoch->sat_[i_sat].sys_ && 
                refsats[i_sys] == sats_epoch->sat_[i_sat].prn_)
                usedsats(i_row, 2) = 1;
        i_row ++;
    }
    return usedsats;
}

void CRtkekf::updateStatus(MatrixXd obs, sat* sats_epoch, int nobs, int* refsats, res_t res, int* sysobs) {
    if (status_ == RTK_INIT) {
        memcpy(refsats_, refsats, sizeof(int) * MAXSYS);
        memcpy(sysobs_, sysobs, sizeof(int) * MAXSYS);
        state_.resize((nobs - opt_->nsys_) * opt_->freqnum_ + 3, 1);
        initvarsys(state_.row(), 1E-8);
        initvarstate(state_.row(), 15 * 15);
        initstate(sats_epoch, res.bpos_ecef_, refsats, res);
        state_trans_.resize(state_.row(), state_.row()); 
        state_trans_.Identity();
        nobs_ = nobs;
        obssats_ = obssats(sats_epoch, nobs, refsats);
        status_ = RTK_NORMAL;
    } else {
        MatrixXd obssat = obssats(sats_epoch, nobs, refsats_);
        if (checkObsChange(obssat, nullptr, nullptr))
            status_ |= RTK_OBSCHANGE;
        for (int i = 0; i < MAXSYS; ++ i) 
            if (refsats_[i] != refsats[i]) status_ |= RTK_REFCHANGE;
        // TODO: WEEK LEAP DETECTION
    }
}

double CRtkekf::dist(double* a, double* b) {
    double dist = 0;
    for (int i = 0; i < 3; ++ i)
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(dist);
}

sat_s CRtkekf::findRef(sat sats, int sysflag, int prn) {
    int nobs = sats.nsats_;
    for (int isat = 0; isat < nobs; ++isat) {
        if (sats.sat_[isat].sys_ == sysflag && sats.sat_[isat].prn_ == prn)
            return sats.sat_[isat];
    }
}

double CRtkekf::GetFreq(int sys, int freqflag) {
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

void CRtkekf::initstate(sat* sats, double* b_pos, int* refsats, res_t res) {
    int ar_pos = 3;
    for (int i = 0; i < 3; ++ i)
        state_(i, 0) = res.rpos_ecef_[i];
    int ref_prn = -1, sysflag = -1;
    int nobs = sats->nsats_;
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
        initambiguity(ref_sat_r, ref_sat_b, sats[1].sat_[isat], sats[0].sat_[isat], ar_pos);
    }
}

void CRtkekf::initambiguity(sat_s ref_sat_r, sat_s ref_sat_b, sat_s sat_r, sat_s sat_b, int &ar_pos) {
    for (int i = 0; i < MAXFREQ; ++i) {
        if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
            double phase_obs = sat_r.obs_->L[i] - ref_sat_r.obs_->L[i]-
                sat_b.obs_->L[i] + ref_sat_b.obs_->L[i];
            double pseu_obs = sat_r.obs_->P[i] - ref_sat_r.obs_->P[i]- 
                sat_b.obs_->P[i] + ref_sat_b.obs_->P[i];
            double freq = GetFreq(sat_r.sys_, FREQ_ARRAY[i]);
            double lambda = VEL_LIGHT / freq;
            state_(ar_pos ++, 0) = (pseu_obs - phase_obs * lambda) / lambda;
        }
    }
}

void CRtkekf::initvarsys(int dim, double var) {
    var_sys_.resize(dim, dim);
    for (int i = 0; i < dim; ++i) 
        for (int j = 0; j < dim; ++ j)
            if (i == j) var_sys_(i, j) = var;
            else        var_sys_(i, j) = 0;
}

void CRtkekf::initvarstate(int dim, double var) {
    var_state_.resize(dim, dim);
    for (int i = 0; i < dim; ++i) 
        for (int j = 0; j < dim; ++ j)
            if (i == j && i < 3) var_state_(i, j) = var;
            else if (i == j && i >= 3) var_state_(i, j) = var / 0.04;
            else        var_state_(i, j) = 0;
}

void CRtkekf::setopt(prcopt* opt) {
    if (!opt) return;
    opt_ = opt;
}

CRtkekf::~CRtkekf() {
    opt_ = nullptr;
}

bool CRtkekf::checkObsChange(MatrixXd obssats_current, MatrixXd* increase, MatrixXd* decrease) {
    // MatrixXd obssat = obssats(sats_epoch, )
    if (checkObsIncrease(obssats_current, increase)) return true;
    if (checkObsDecrease(obssats_current, decrease)) return true;

    return false;
}

bool CRtkekf::checkObsIncrease(MatrixXd obssats_current, MatrixXd* increase) {
    if (increase) {
        increase->resize(MAX_OBSCHANGE, MAXSYS);
        increase->Zero();
    }
    int sys_change[MAXSYS] = {0};
    bool isincrease = false;
    for (int i = 0; i < obssats_current.row(); ++ i) {
        bool isfind = false;
        for (int j = 0; j < obssats_.row(); ++ j) {
            if (obssats_current(i, 0) == obssats_(j, 0) && 
                obssats_current(i, 1) == obssats_(j, 1)) {
                isfind = true; break;
            }
        }
        if (isfind) continue;
        isincrease = true;
        if (!increase) continue; 
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys)
            if (SYS_ARRAY[i_sys] == int(obssats_current(i, 0)))
                (*increase)(sys_change[i_sys] ++, i_sys) = obssats_current(i, 1);
    }
    return isincrease;
}

bool CRtkekf::checkObsDecrease(MatrixXd obssats_current, MatrixXd* decrease) {
    if (decrease) {
        decrease->resize(MAX_OBSCHANGE, MAXSYS);
        decrease->Zero();
    }
    int sys_change[MAXSYS] = {0};
    bool isdecrease = false;
    for (int i = 0; i < obssats_.row(); ++ i) {
        bool isfind = false;
        for (int j = 0; j < obssats_current.row(); ++ j) {
            if (obssats_current(j, 0) == obssats_(i, 0) && 
                obssats_current(j, 1) == obssats_(i, 1)) {
                isfind = true; break;
            }
        }
        if (isfind) continue;
        isdecrease = true;
        if (!decrease) continue; 
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys)
            if (SYS_ARRAY[i_sys] == int(obssats_(i, 0)))
                (*decrease)(sys_change[i_sys] ++, i_sys) = obssats_(i, 1);
    }
    return isdecrease;
}


void CRtkekf::statusFix(sat* sats_epoch, int nobs, int* refsats, int* obssys) {
    for (int i = 0; i < RTK_STATUS_NUM; ++i) {
        if ((status_ & STATUS_ARRAY[i]) == STATUS_ARRAY[i]) 
            switch (STATUS_ARRAY[i]) {
            case RTK_LEAP:
                break;
            case RTK_REFCHANGE:
                refchange(refsats, obssys); status_ &= 0x3B;
                memcpy(refsats_, refsats, sizeof(int) * MAXSYS);
                // obssats_ = obssats(sats_epoch, nobs, refsats_);
                break;
            case RTK_OBSCHANGE: {
                MatrixXd obssats_epoch = obssats(sats_epoch, nobs, refsats_);
                MatrixXd increase;
                MatrixXd decrease;
                if (checkObsIncrease(obssats_epoch, &increase)) {
                    fixNumIncrease(sats_epoch, obssys, nobs, increase);
                    // status_ = RTK_NORMAL;
                    nobs_ = nobs;
                }
                if (checkObsDecrease(obssats_epoch, &decrease)) {
                    fixNumDecrease(sats_epoch, obssys, nobs, decrease);
                    // status_ = RTK_NORMAL;
                    nobs_ = nobs;
                }
                status_ = RTK_NORMAL;
                cout << increase << endl << endl;
                cout << decrease << endl << endl;
                break;
            }
            default:
                break;
            }
    }
    memcpy(sysobs_, obssys, sizeof(int) * MAXSYS);
    // ofstream out("./kf_debug.txt");
    // out << "*****************************\n" << obssats_ << endl << endl;
    obssats_ = obssats(sats_epoch, nobs, refsats_);
    // out << obssats_ << endl << endl;
    // out.close();

}

void CRtkekf::fixNumIncrease(sat* sats, int* obssys, int nobs, MatrixXd increase) {
    int pos = -1, num_sys_obs = 0;
    int freqnum = opt_->freqnum_;
    bool isAfterRef = false;
    ofstream out ("./kf_debug.txt");
    out << increase << endl << endl;
    out << obssats_ << endl << endl;
    auto test = obssats(sats, nobs, refsats_);
    out << test << endl << endl;
    for (int i = 0; i < sats->nsats_; ++ i) {
        if (!sats->sat_[i].isused) continue;
        bool isref = false;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (SYS_ARRAY[i_sys] == sats[1].sat_[i].sys_ &&
                sats[1].sat_[i].prn_ == refsats_[i_sys])
                isref = true;
        }
        if (isref) {
            isAfterRef = true;
            continue;
        } 
        pos += 1;
        if (!havesat(sats->sat_[i], increase)) continue;

        // construct transfer matrix
        MatrixXd T(var_state_.row() + freqnum, var_state_.col());
        MatrixXd T_info(obssats_.row() + 1, obssats_.row());  // transfer matrix info
        T.Zero(); T_info.Zero();
        for (int j = 0; j < 3; ++ j) T(j, j) = 1;
        int row = 3, col = 3;
        int sysflag = 0;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys)
            if (sats->sat_[i].sys_ == SYS_ARRAY[i_sys])
                sysflag = i_sys;
        int total_num = 0;
        for (int tmp = 0; tmp < sysflag; ++ tmp) {
            total_num += (tmp == sysflag ? sysobs_[tmp] : sysobs_[tmp] - 1); 
        }
        int pos_add = pos - total_num;
        for (int i_part = 0; i_part < opt_->nsys_ * opt_->freqnum_;  ++ i_part) {

            // for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                // if (sats->sat_[i].sys_ != SYS_ARRAY[i_sys]) continue;
            int sys_dd_num = sysobs_[int(i_part / freqnum)]; // contain new sats and ref sats
            if (i_part == sysflag * 2 || i_part - 1 == sysflag * 2) {
                // if (i_part - 1 == sysflag * 2) pos_add -= 1;  // new line added at i_part == i_sys * 2
                MatrixXd T_part(sys_dd_num, sys_dd_num - 1); T_part.Zero();
                for (int i_row_t = 0; i_row_t < T_part.row(); ++ i_row_t) {
                    // int minus = i_part - sysflag * 2;
                    if (i_row_t< pos_add) T_part(i_row_t, i_row_t) = 1;
                    else if (i_row_t == pos_add) continue;
                    else    T_part(i_row_t, i_row_t - 1) = 1;
                }
                out << T_part << endl << endl;
                T.block(row, col, T_part);
                row += T_part.row(); col += T_part.col();
            } else if (i_part != sysflag * 2 || i_part - 1 != sysflag * 2) {
                MatrixXd T_part(sys_dd_num - 1, sys_dd_num - 1);    // minus 1 to exclude ref sat
                T_part.Identity();
                T.block(row, col, T_part);
                row += T_part.row(); col += T_part.col();
            }

            // }
        }
        sysobs_[sysflag] += 1;
        out << T << endl << endl;
        out << state_ << endl << endl;
        state_ = T * state_;
        int ddobs[MAXSYS] = {0};
        initstate_s(sats, sats->sat_[i].prn_, sats[1].sat_[i].sys_, ddobs);
        for (int j = 0, num = 0; j < state_.row(); ++ j)
            if (state_(j, 0) == 0)
                state_(j, 0) = ddobs[num ++];
        out << state_ << endl << endl;
        var_state_ = T * var_state_ * T.transpose();
        state_trans_.resize(state_.row(), state_.row()); state_trans_.Identity();
        var_sys_.resize(state_.row(), state_.row()); var_sys_.Zero();

        if (isAfterRef) pos += 1; // no need to consider ref sats for obssats
        pos += sysflag;
        for (int i_row = 0; i_row < T_info.row(); ++i_row) {
            if (i_row == pos) continue;
            if (i_row < pos)
                T_info (i_row, i_row) = 1;
            else 
                T_info (i_row, i_row - 1) = 1;   // dual-freq
        }
        out << obssats_ << endl;
        auto test = obssats(sats, nobs, refsats_);
        out << "test\n" << test << endl << endl;
        obssats_ = T_info * obssats_;
        for (int i_row = 0; i_row < obssats_.row(); ++ i_row)
            if (obssats_(i_row, 0) == 0) {
                obssats_(i_row, 0) = sats->sat_[i].sys_;
                obssats_(i_row, 1) = sats->sat_[i].prn_;
            }
        for (int i_row = 0; i_row < var_state_.row(); ++i_row) {
            if (abs(var_state_(i_row, i_row)) <= 1e-8)
                var_state_(i_row, i_row) = INIT_STATE_VAR * INIT_STATE_VAR;
        }
        out << var_state_ << endl << endl;
        // for (int i_row = 0; i_row < state_.row(); ++ i_row) 
        //     if (state_(i_row, 0) == 0) 
        //         initstate_s(sats, i_row);
        out << state_ << endl << endl;
        out << var_sys_ << endl;
        out << obssats_ << endl << endl;
    }
    nobs_ = nobs;
    out.close();
}

void CRtkekf::fixNumDecrease(sat* sats, int* obssys, int nobs, MatrixXd decrease) {
    int num_sys_obs = 0;
    int freqnum = opt_->freqnum_;
    bool isafterRef = false;
    ofstream out ("./kf_debug.txt");
    out << decrease << endl << endl;
    out << obssats_ << endl << endl;
    int ref_pos[MAXSYS];
    for (int i = 0; i < MAXSYS; i ++) ref_pos[i] = -1;
    for (int i = 0; i < obssats_.row(); ++ i) {
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (refsats_[i_sys] == obssats_(i, 1) && 
                SYS_ARRAY[i_sys] == obssats_(i, 0) && 
                obssats_(i, 2) == 1)
                ref_pos[i_sys] = i;
        }
    }

    for (int i_row_c = 0; i_row_c < decrease.row(); ++ i_row_c) {
        for (int i_row_o = 0; i_row_o < obssats_.row(); ++ i_row_o) {
            for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                if (SYS_ARRAY[i_sys] != obssats_(i_row_o, 0) ||
                    obssats_(i_row_o, 1) != decrease(i_row_c, i_sys))
                    continue;
                MatrixXd T_trans(state_.row() - freqnum, state_.row()); T_trans.Zero();
                for (int i = 0; i < 3; ++ i) T_trans(i, i) = 1;
                int pos = i_row_o, total_obs = 0;
                int row = 3, col = 3;
                for (int i = 0; i < i_sys; ++ i) total_obs += sysobs_[i] - 1;
                pos -= (total_obs == 0 ? total_obs : total_obs + 1);
                ref_pos[i_sys] -= (total_obs == 0 ? total_obs : total_obs + 1);
                if (pos >= ref_pos[i_sys]) {
                    pos -= 1; isafterRef = true;
                }
                for (int i_part = 0; i_part < opt_->freqnum_ * opt_->nsys_; ++ i_part) {
                    if (i_part == i_sys * 2 || i_part - 1 == i_sys * 2) {
                        int ddobs = sysobs_[i_part / freqnum] - 1;
                        MatrixXd T_part (ddobs - 1, ddobs); T_part.Zero();
                        for (int i = 0; i < T_part.row(); ++ i) 
                            if (i >= pos) T_part(i, i + 1) = 1;
                            else T_part(i, i) = 1;
                        T_trans.block(row, col, T_part);
                        row += T_part.row(); col += T_part.col();
                    } else if (i_part != i_sys * 2 || i_part - 1 != i_sys * 2) {
                        int ddobs = sysobs_[i_part / freqnum] - 1;
                        MatrixXd T_part (ddobs, ddobs); T_part.Identity();
                        T_trans.block(row, col, T_part);
                        row += T_part.row(); col += T_part.col();
                    }

                }
                sysobs_[i_sys] -= 1;
                out << state_ << endl;
                out << T_trans << endl << endl;
                state_ = T_trans * state_;
                var_state_ = T_trans * var_state_ * T_trans.transpose();
                out << var_state_ << endl;
                if (isafterRef) pos += 1;
                pos += (total_obs == 0 ? total_obs : total_obs + 1);
                MatrixXd T_info(obssats_.row() - 1, obssats_.row()); T_info.Zero();
                for (int i = 0; i < T_info.row(); ++ i) {
                    if (i < pos) T_info(i, i) = 1;
                    else T_info(i, i + 1) = 1;
                }
                obssats_ = T_info * obssats_;
                state_trans_.resize(state_.row(), state_.row()); state_trans_.Identity();
                var_sys_.resize(state_.row(), state_.row()); var_sys_.Zero();
                nobs_ = nobs;
                out << obssats_ << endl << endl;
                out << state_ << endl << endl;
                out << T_info << endl << endl;
                out << T_trans << endl << endl;  
            }
        }
    }
}


bool CRtkekf::havesat(sat_s sat, MatrixXd changeObs) {
    for (int i_row = 0; i_row < MAX_OBSCHANGE; ++ i_row) {
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (SYS_ARRAY[i_sys] == sat.sys_ && 
                sat.prn_ == changeObs(i_row, i_sys)) 
                return true;
        }
    }
    return false;
}

bool CRtkekf::havesat(MatrixXd obssats, MatrixXd changeObs, int &pos) {
    pos = 0;
    for (int i_row_c = 0; i_row_c < changeObs.row(); ++ i_row_c) {
        for (int i_row_o = 0; i_row_o < obssats.row(); ++ i_row_o) {
            pos += 1;
            for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                if (SYS_ARRAY[i_sys] == obssats(i_row_o, 0) &&
                    obssats(i_row_o, 1) == changeObs(i_row_c, i_sys))
                    return true;
            }
        }
    }
    return false;
}


void CRtkekf::refchange(int* refsats, int* obssys) {
    int new_ref_position[MAXSYS];
    for (int i = 0; i < MAXSYS; ++ i) new_ref_position[i] = -1;
    int freqnum = opt_->freqnum_;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i]);

        for (int i_row = 0; i_row < obssats_.row(); ++ i_row) 
            if (obssats_(i_row, 0) == SYS_ARRAY[i_sys] && 
                obssats_(i_row, 1) == refsats[i_sys] && 
                obssats_(i_row, 2) != 1) {
                // position of L1-ambiguity 
                new_ref_position[i_sys] = i_row - total; 
            }
        
        for (int i_row = 0; i_row < obssats_.row(); ++ i_row) 
            if (obssats_(i_row, 0) == SYS_ARRAY[i_sys] && 
                obssats_(i_row, 2) == 1)
                if (i_row <= new_ref_position[i_sys] + total)
                    new_ref_position[i_sys] -= 1;
        // if (total != 0) new_ref_position[i_sys] -= 1;
    }
    
    MatrixXd trans(state_.row(), state_.row()); trans.Zero();
    for (int i = 0; i < 3; ++ i) trans(i, i) = 1;

    int row = 3;
    // for (int i_part = 0; i_part < opt_->freqnum_ * opt_->nsys_; ++ i_part)
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i] - 1) * freqnum;
        int ddobs = obssys[i_sys] - 1;
        if (refsats_[i_sys] == refsats[i_sys] || new_ref_position[i_sys] == -1) {
            MatrixXd identity(ddobs * freqnum, ddobs * freqnum);
            identity.Identity();
            for (int i_freq = 0; i_freq < freqnum; ++ i_freq)
                trans.block(row + total + i_freq * obssys[i_sys], row + total + i_freq * obssys[i_sys], identity);
            continue;
        }
        MatrixXd trans_part(ddobs, ddobs);

        for (int ifreq = 0; ifreq < opt_->freqnum_; ++ ifreq) {
            trans_part.Zero();
            for (int i = 0; i < trans_part.row(); ++ i) {
                trans_part(i, new_ref_position[i_sys]) = -1;
                if (i == new_ref_position[i_sys]) continue;
                trans_part(i, i) = 1;
            }
            cout << trans_part << endl << endl;
            trans.block(row + total + ddobs * ifreq, row + total + ddobs * ifreq, trans_part);
        }
        
        // for (int i = i_sys, total = 0; i >= 0; i --) total += (obssys[i] - 1) * freqnum;
        // for (; row < trans.row() && row < total + 3; row += total) {
        //     for (int i_freq = 0; i_freq < freqnum; ++ i_freq) {
        //         trans(row, new_ref_position[i_sys] + i_freq * total / freqnum) = -1;
        //         if (row == new_ref_position[i_sys] + i_freq * total / freqnum) continue;
        //         trans(row + i_freq * total / freqnum, row + i_freq * total / freqnum) = 1; // TODO: not tested
        //     }
        // }
    }
    ofstream out ("./trans.txt");
    out << obssats_ << endl << endl;
    out << refsats[0] << "  " << refsats[1] << endl << endl;
    out << state_ << endl << endl;
    state_ = trans * state_;
    var_state_ = trans * var_state_ * trans.transpose();
    out << state_ << endl << endl;
    out << trans << endl << endl;
    out.close();
}

void CRtkekf::initstate_s(sat* sats, int prn, int sys, int* ddobs) {
    int nobs = sats->nsats_;
    int ref_prn = -1, sysflag = -1;
    int freqnum = opt_->freqnum_;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = state_(i, 0);
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((opt_->navsys_ & SYS_ARRAY[i]) == sats[1].sat_[isat].sys_) {
                ref_prn = refsats_[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if (sats->sat_[isat].sys_ != sys || sats->sat_[isat].prn_ != prn)
            continue;
        // if(sats[1].sat_[isat].prn_ == ref_prn && sats[1].sat_[isat].sys_ == sysflag) {
        //     pos += 2;
        //     continue;
        // }
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        sat_s sat_r = sats[1].sat_[isat], sat_b = sats[0].sat_[isat];
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sat_r.obs_->L[i] - ref_sat_r.obs_->L[i]-
                    sat_b.obs_->L[i] + ref_sat_b.obs_->L[i];
                double pseu_obs = sat_r.obs_->P[i] - ref_sat_r.obs_->P[i]- 
                    sat_b.obs_->P[i] + ref_sat_b.obs_->P[i];
                double freq = GetFreq(sat_r.sys_, FREQ_ARRAY[i]);
                double lambda = VEL_LIGHT / freq;
                ddobs[num ++] = (pseu_obs - phase_obs * lambda) / lambda;
            }
        }
        return;
    }
}

