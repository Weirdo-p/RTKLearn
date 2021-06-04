#include "navigation/optimal/rtkekf.h"
#include "navigation/pnt/pntbase.h"

bool CRtkekf::optimize(sat* sats_epoch, res_t &res) {
    int refsats[MAXSYS];
    int nobs = obsnumber(sats_epoch);
    chooseref(sats_epoch[1], refsats);
    int row = 0, col = 0;
    getDesignDim(sats_epoch[1], nobs, row, col);
    int sysobs[MAXSYS] = {0};
    getSysObs(sats_epoch[1], sysobs);
    design_.resize(row, col); design_.Zero();
    var_obs_.resize(row, row);
    MatrixXd obs(row, 1);
    getDesign(sats_epoch, nobs, res.rpos_ecef_, refsats, design_);
    getweight(sats_epoch, refsats, nobs, var_obs_);
    getddobs(sats_epoch, res.rpos_ecef_, refsats, sysobs, res, obs);
    int issuccess = false;
    var_obs_ = var_obs_.inverse(issuccess);
    updateStatus(obs, sats_epoch, nobs, refsats, res, sysobs);
    statusFix(sats_epoch, nobs, refsats, sysobs);
    gfcycle(sats_epoch);
    issuccess = 0;
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
    return true;
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
                double freq = CPntbase::GetFreq(sysflag, FREQ_ARRAY[i_freq]);
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

CRtkekf::CRtkekf(prcopt* opt) {
    if (!opt_) {
        opt_ = new prcopt;
        memcpy(opt_, opt, sizeof(prcopt));
        opt_->nbase_ = opt->nbase_; opt_->nrover_ = opt->nrover_;
    }
    memset(LeapSatsPos_, 0, sizeof(int) * MAXOBS);
    status_ = RTK_INIT;
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

void CRtkekf::initstate(sat* sats, double* b_pos, int* refsats, res_t res) {
    int ar_pos = 3;
    for (int i = 0; i < 3; ++ i)
        state_(i, 0) = res.rpos_ecef_[i];
    int ref_prn = -1, sysflag = -1, sys_pos = -1;
    int nobs = sats->nsats_;
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1].sat_[isat].isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((opt_->navsys_ & SYS_ARRAY[i]) == sats[1].sat_[isat].sys_) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i]; sys_pos = i;
                break;
            }
        }
        if(sats[1].sat_[isat].prn_ == ref_prn && sats[1].sat_[isat].sys_ == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        int nddobs = sysobs_[sys_pos] - 1;
        initambiguity(ref_sat_r, ref_sat_b, sats[1].sat_[isat], sats[0].sat_[isat], nddobs, ar_pos);
    }
}

void CRtkekf::initambiguity(sat_s ref_sat_r, sat_s ref_sat_b, sat_s sat_r, sat_s sat_b, int nddobs, int &ar_pos) {
    int num = 0;
    for (int i = 0; i < MAXFREQ; ++i) {
        if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
            double phase_obs = sat_r.obs_->L[i] - ref_sat_r.obs_->L[i]-
                sat_b.obs_->L[i] + ref_sat_b.obs_->L[i];
            double pseu_obs = sat_r.obs_->P[i] - ref_sat_r.obs_->P[i]- 
                sat_b.obs_->P[i] + ref_sat_b.obs_->P[i];
            double freq = CPntbase::GetFreq(sat_r.sys_, FREQ_ARRAY[i]);
            double lambda = VEL_LIGHT / freq;
            if (num == 0)
                state_(ar_pos, 0) = -(pseu_obs - phase_obs * lambda) / lambda;
            else 
                state_(ar_pos + nddobs, 0) = -(pseu_obs - phase_obs * lambda) / lambda;
            ++ num;
        }
    }
    ar_pos += 1;
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

                break;
            case RTK_OBSCHANGE: {
                MatrixXd obssats_epoch = obssats(sats_epoch, nobs, refsats_);
                MatrixXd increase;
                MatrixXd decrease;
                if (checkObsIncrease(obssats_epoch, &increase)) {
                    fixNumIncrease(sats_epoch, obssys, nobs, increase);
                    nobs_ = nobs;
                }
                if (checkObsDecrease(obssats_epoch, &decrease)) {
                    fixNumDecrease(sats_epoch, obssys, nobs, decrease);
                    nobs_ = nobs;
                }
                status_ = RTK_NORMAL;
                break;
            }
            default:
                break;
            }
    }
    memcpy(sysobs_, obssys, sizeof(int) * MAXSYS);
    obssats_ = obssats(sats_epoch, nobs, refsats_);
    status_ = RTK_NORMAL;
}

void CRtkekf::fixNumIncrease(sat* sats, int* obssys, int nobs, MatrixXd increase) {
    int pos = -1, num_sys_obs = 0;
    int freqnum = opt_->freqnum_;
    bool isAfterRef = false;

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
            total_num += sysobs_[tmp] - 1; 
        }
        int pos_add = pos - total_num;
        // if (isAfterRef) pos_add -= 1;
        for (int i_part = 0; i_part < opt_->nsys_ * opt_->freqnum_;  ++ i_part) {
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

                T.block(row, col, T_part);
                row += T_part.row(); col += T_part.col();
            } else if (i_part != sysflag * 2 || i_part - 1 != sysflag * 2) {
                MatrixXd T_part(sys_dd_num - 1, sys_dd_num - 1);    // minus 1 to exclude ref sat
                T_part.Identity();
                T.block(row, col, T_part);
                row += T_part.row(); col += T_part.col();
            }
        }
        sysobs_[sysflag] += 1;
        state_ = T * state_;
        double ddobs[MAXSYS] = {0};
        ofstream out("./debug_increase.txt");
        out << obssats_ << endl << endl;
        out << increase << endl << endl;
        initstate_s(sats, sats->sat_[i].prn_, sats[1].sat_[i].sys_, ddobs);
        out << state_ << endl << endl;
        out << ddobs[0] << "  " << ddobs[1] << endl << endl;
        for (int j = 0, num = 0; j < state_.row(); ++ j)
            if (state_(j, 0) == 0)
                state_(j, 0) = ddobs[num ++];

        var_state_ = T * var_state_ * T.transpose();
        state_trans_.resize(state_.row(), state_.row()); state_trans_.Identity();
        var_sys_.resize(state_.row(), state_.row()); var_sys_.Zero();
        out << T << endl << endl;
        out << state_ << endl << endl;
        out << var_state_ << endl << endl;
        out.close();
        int tmp = pos;
        if (isAfterRef) pos += 1; // no need to consider ref sats for obssats
        pos += sysflag;
        for (int i_row = 0; i_row < T_info.row(); ++i_row) {
            if (i_row == pos) continue;
            if (i_row < pos)
                T_info (i_row, i_row) = 1;
            else 
                T_info (i_row, i_row - 1) = 1;   // dual-freq
        }
        obssats_ = T_info * obssats_;
        for (int i_row = 0; i_row < obssats_.row(); ++ i_row)
            if (obssats_(i_row, 0) == 0) {
                obssats_(i_row, 0) = sats->sat_[i].sys_;
                obssats_(i_row, 1) = sats->sat_[i].prn_;
            }
        for (int i_row = 0; i_row < var_state_.row(); ++i_row) {
            if (abs(var_state_(i_row, i_row)) <= 1e-8)
                var_state_(i_row, i_row) = INIT_STATE_VAR * INIT_STATE_VAR / 0.19 * 100;
        }
        pos = tmp;  // for any number of satellites increase
    }
    nobs_ = nobs;
}

void CRtkekf::fixNumDecrease(sat* sats, int* obssys, int nobs, MatrixXd decrease) {
    int num_sys_obs = 0;
    int freqnum = opt_->freqnum_;
    bool isafterRef = false;
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
    ofstream out ("./debug_decrease.txt");
    out << decrease << endl << endl;
    out << obssats_ << endl << endl;
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
                out << T_trans << endl << endl;
                sysobs_[i_sys] -= 1;
                out << state_ << endl << endl;
                state_ = T_trans * state_;
                var_state_ = T_trans * var_state_ * T_trans.transpose();
                out << state_ << endl << endl;
                if (isafterRef) pos += 1;
                pos += (total_obs == 0 ? total_obs : total_obs + 1);
                MatrixXd T_info(obssats_.row() - 1, obssats_.row()); T_info.Zero();
                for (int i = 0; i < T_info.row(); ++ i) {
                    if (i < pos) T_info(i, i) = 1;
                    else T_info(i, i + 1) = 1;
                }
                obssats_ = T_info * obssats_;
                out << obssats_ << endl << endl;
                state_trans_.resize(state_.row(), state_.row()); state_trans_.Identity();
                var_sys_.resize(state_.row(), state_.row()); var_sys_.Zero();
                nobs_ = nobs;
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
    }
    
    MatrixXd trans(state_.row(), state_.row()); trans.Zero();
    for (int i = 0; i < 3; ++ i) trans(i, i) = 1;

    int row = 3;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i] - 1) * freqnum;
        int ddobs = obssys[i_sys] - 1;
        if (refsats_[i_sys] == refsats[i_sys] || new_ref_position[i_sys] == -1) {
            MatrixXd identity(ddobs, ddobs);
            identity.Identity();
            for (int i_freq = 0; i_freq < freqnum; ++ i_freq)
                trans.block(row + total + i_freq * (obssys[i_sys] - 1), row + total + i_freq * (obssys[i_sys] - 1), identity);
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
            trans.block(row + total + ddobs * ifreq, row + total + ddobs * ifreq, trans_part);
        }
        
    }
    state_ = trans * state_;
    var_state_ = trans * var_state_ * trans.transpose();
    sortRefChange(new_ref_position, refsats, obssys);
}

void CRtkekf::sortRefChange(int* new_ref_position, int* refsats, int* obssys) {
    int old_ref_position[MAXSYS];
    int freqnum = opt_->freqnum_;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        if (new_ref_position[i_sys] == -1) continue;
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += obssys[i];

        for (int i_row = 0; i_row < obssats_.row(); ++ i_row) 
            if (obssats_(i_row, 0) == SYS_ARRAY[i_sys] && 
                obssats_(i_row, 1) == refsats_[i_sys]) {
                // position of L1-ambiguity 
                old_ref_position[i_sys] = i_row - total; 
            }
        if (old_ref_position[i_sys] > new_ref_position[i_sys])
            old_ref_position[i_sys] -= 1;
        
        total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i] - 1);

        int ddobs = obssys[i_sys] - 1;
        for (int i_freq = 0; i_freq < freqnum; ++ i_freq) {
            int row = 3 + total * freqnum + ddobs * i_freq;
            // for the part will be sorted to obssats_ position
            MatrixXd part_state(ddobs, state_.col());
            MatrixXd part_var_state(ddobs, ddobs);

            for (int i = 0; i < ddobs; ++ i)
                part_state(i, 0) = state_(row + i, 0);
            for (int i = 0; i < ddobs; ++ i)
                for (int j = 0; j < part_var_state.col(); ++ j)
                    part_var_state(i, j) = var_state_(row + i, row + j);
            
            MatrixXd temp_state(1, state_.col());
            int new_ref_pos = new_ref_position[i_sys];
            int old_ref_pos = old_ref_position[i_sys];
            temp_state(0, 0) = part_state(new_ref_pos, 0);
            MatrixXd T_part(part_state.row(), part_state.row()); T_part.Zero();
            MatrixXd T(var_state_.row(), var_state_.col());

            T_part(old_ref_pos, new_ref_pos) = 1;

            if (new_ref_pos == old_ref_pos) continue;
            else if (new_ref_pos < old_ref_pos) // up move
                for (int i = 0; i < old_ref_pos - new_ref_pos; ++ i) {
                    part_state(new_ref_pos + i, 0) = part_state(new_ref_pos + i + 1, 0);
                    T_part(new_ref_pos + i, new_ref_pos + i + 1) = 1;
                }
            else // down move
                for (int i = new_ref_pos - 1; i >= old_ref_pos; -- i) {
                    part_state(i + 1, 0) = part_state(i, 0);
                    T_part(i + 1, i) = 1;
                }
            part_state(old_ref_pos, 0) = temp_state(0, 0);

            for (int i = 0; i < T_part.row(); ++ i) {
                bool iszero = true;
                for (int j = 0; j < T_part.col(); ++ j)
                    if (T_part(i, j) != 0) iszero = false;
                if (!iszero) continue;
                T_part(i, i) = 1;
            }
            T.Identity(); T.block(row, row, T_part);
            var_state_ = T * var_state_ * T.transpose();
            state_.block(row, 0, part_state);
        }
    }

}

void CRtkekf::initstate_s(sat* sats, int prn, int sys, double* ddobs) {
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
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        sat_s sat_r = sats[1].sat_[isat], sat_b = sats[0].sat_[isat];
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sat_r.obs_->L[i] - ref_sat_r.obs_->L[i]-
                    sat_b.obs_->L[i] + ref_sat_b.obs_->L[i];
                double pseu_obs = sat_r.obs_->P[i] - ref_sat_r.obs_->P[i]- 
                    sat_b.obs_->P[i] + ref_sat_b.obs_->P[i];
                double freq = CPntbase::GetFreq(sat_r.sys_, FREQ_ARRAY[i]);
                double lambda = VEL_LIGHT / freq;
                ddobs[num ++] = -(pseu_obs - phase_obs * lambda) / lambda;
            }
        }
        return;
    }
}

void CRtkekf::gfcycle(sat* sats_epoch) {
    static sat* sat_last_epoch = nullptr;
    const double THRESHOLD = 0.0131; // unit: meter
    if (opt_->freqnum_ < 2) return; // at least dual-freq
    for (int i_site = 0; i_site < opt_->sitenum_; ++ i_site) {
        int pos = -1;
        bool isafterref = false;
        for (int i_sat = 0; i_sat < sats_epoch[i_site].nsats_; ++ i_sat) {
            if (!sats_epoch[i_site].sat_[i_sat].isused) continue;
            pos += 1;

            // gf combination
            double gf = 0;
            for (int i_freq = 0; i_freq < MAXFREQ; ++ i_freq) {
                if ((opt_->freqtype_ & FREQ_ARRAY[i_freq]) == FREQ_ARRAY[i_freq]) {
                    double freq = CPntbase::GetFreq(sats_epoch[i_site].sat_[i_sat].sys_, FREQ_ARRAY[i_freq]);
                    double lambda = VEL_LIGHT / freq;
                    gf += sats_epoch[i_site].sat_[i_sat].obs_->L[i_freq] * lambda * pow(-1, i_freq);
                }
            }
            sats_epoch[i_site].sat_[i_sat].gf_ = gf;
            if (!sat_last_epoch) continue; // not initialized 
            int sysflag = 0;
            for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) 
                if (SYS_ARRAY[i_sys] == sats_epoch[i_site].sat_[i_sat].sys_)
                    sysflag = i_sys;
            if (sats_epoch[i_site].sat_[i_sat].prn_ == refsats_[sysflag]) isafterref = true;
            int state_pos = pos;
            int ddobs = 0, freqnum = opt_->freqnum_;
            for (int i_sys = 0; i_sys < sysflag; ++ i_sys) ddobs += (sysobs_[i_sys] - 1);
            state_pos -= (ddobs + sysflag); state_pos += freqnum * ddobs;
            if (isafterref) state_pos -= 1;
            state_pos += 3;
            double last_gf = searchgf(sat_last_epoch[i_site], sats_epoch[i_site].sat_[i_sat]);
            if (last_gf == -1) continue;

            if (abs(gf - last_gf) >= THRESHOLD) { // cycle slip occurred, reset state
                auto ref_sat_b = findRef(sats_epoch[0], sats_epoch[0].sat_[i_sat].sys_, refsats_[sysflag]);
                auto ref_sat_r = findRef(sats_epoch[1], sats_epoch[1].sat_[i_sat].sys_, refsats_[sysflag]);
                auto sat_b = sats_epoch[0].sat_[i_sat];
                auto sat_r = sats_epoch[1].sat_[i_sat];
                initambiguity(ref_sat_r, ref_sat_b, sat_r, sat_b, sysobs_[sysflag] - 1, state_pos);
                state_pos -= 1;
                for (int i = 0; i < var_state_.row(); ++ i) 
                    if (i != state_pos) var_state_(state_pos, i) = 0;
                    else var_state_(state_pos, state_pos) = 30 * 30;
                state_pos += sysobs_[sysflag] - 1;

                for (int i = 0; i < var_state_.row(); ++ i) 
                    if (i != state_pos) var_state_(state_pos, i) = 0;
                    else var_state_(state_pos, state_pos) = 30 * 30;
            }
        } // each satellite loop
    } // each site loop

    if (!sat_last_epoch) sat_last_epoch = new sat[opt_->sitenum_];
    memcpy(sat_last_epoch, sats_epoch, sizeof(sat) * opt_->sitenum_);
}

double CRtkekf::searchgf(sat sat_last, sat_s obj) {
    for (int isat = 0; isat < sat_last.nsats_; ++ isat)
        if (sat_last.sat_[isat].sys_ == obj.sys_ && sat_last.sat_[isat].prn_ == obj.prn_)
            return sat_last.sat_[isat].gf_;
    return -1;
}

int CRtkekf::obsnumber(sat* sats) {
    int nsites = opt_->sitenum_;
    int nobs = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < sats->nsats_; ++ i_sat) {
        if (!sats[0].sat_[i_sat].isused || !sats[1].sat_[i_sat].isused)
            sats[0].sat_[i_sat].isused = sats[1].sat_[i_sat].isused = false;
        else nobs += 1;
    }
    return nobs;
}

void CRtkekf::getDesignDim(sat sats, int nobs, int &row, int &col) {
    int nsys = opt_->nsys_, nfreq = opt_->freqnum_;
    row = col = 0;    
    row = (nobs - opt_->nsys_) * nfreq * 2;
    col = (nobs - opt_->nsys_) * nfreq + 3;
}

void CRtkekf::getddobs(sat* sats, double* sitepos, int* refsats, int* sysobs, res_t res, MatrixXd &w) {
    int nobs = sats->nsats_;
    int ref_prn = -1, sysflag = -1;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = res.rpos_ecef_[i];
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
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((opt_->freqtype_ & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sats[1].sat_[isat].obs_->L[i] - ref_sat_r.obs_->L[i]-
                    sats[0].sat_[isat].obs_->L[i] + ref_sat_b.obs_->L[i];
                double pseu_obs = sats[1].sat_[isat].obs_->P[i] - ref_sat_r.obs_->P[i]- 
                    sats[0].sat_[isat].obs_->P[i] + ref_sat_b.obs_->P[i];
                double freq = CPntbase::GetFreq(sysflag, FREQ_ARRAY[i]);
                w(num ++, 0) = pseu_obs;
                w(num ++, 0) = phase_obs * (VEL_LIGHT / freq);
            }
        }
    }
}

MatrixXd CRtkekf::GetState() {
    return state_;
}

MatrixXd CRtkekf::GetVar() {
    return var_state_;
}

double CRtkekf::GetInternalSigma() {
    return 1.0;
}
