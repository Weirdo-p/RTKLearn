#include "navigation/optimal/rtkekf.h"
#include "navigation/pnt/pntbase.h"

bool CRtkekf::optimize(sat* sats_epoch, res_t &res) {
    // ofstream out("debug.txt");
    // out << _obssats << endl << endl;
    int refsats[MAXSYS];
    int nobs = obsnumber(sats_epoch);
    chooseref(sats_epoch[1], refsats);
    int row = 0, col = 0;
    getDesignDim(sats_epoch[1], nobs, row, col);
    int sysobs[MAXSYS] = {0};
    getSysObs(sats_epoch[1], sysobs);
    _design.resize(row, col); _design.Zero();
    _var_obs.resize(row, row);
    MatrixXd obs(row, 1);
    getDesign(sats_epoch, nobs, res._rpos_ecef, refsats, _design);
    getweight(sats_epoch, refsats, nobs, _var_obs);
    getddobs(sats_epoch, res._rpos_ecef, refsats, sysobs, res, obs);
    int issuccess = false;
    _var_obs = _var_obs.inverse(issuccess);
    updateStatus(obs, sats_epoch, nobs, refsats, res, sysobs);
    statusFix(sats_epoch, nobs, refsats, sysobs);
    gfcycle(sats_epoch);
    issuccess = 0;
    MatrixXd state_time_predict = _state_trans * _state;
    MatrixXd var_state_time_predict = _state_trans * _var_state * _state_trans.transpose() + _var_sys;
    MatrixXd approx(obs.row(), 1);
    getDDApprox(sats_epoch, state_time_predict, res._bpos_ecef, refsats, approx);
    MatrixXd K = var_state_time_predict * _design.transpose() * (
        _design * var_state_time_predict * _design.transpose() + _var_obs
    ).inverse(issuccess);
    _state = state_time_predict + K * (obs - approx);
    // auto tmp = obs - approx;
    // out << _obssats << endl << endl;
    // cout << tmp << endl << endl;
    
    MatrixXd temp = K * _design;
    MatrixXd identity(temp.row(), temp.col()); identity.Identity();
    _var_state = (identity - temp) * var_state_time_predict * (identity - temp).transpose() + K * _var_obs * K.transpose();
    for (int i = 0; i < 3; ++ i)
        res._rpos_ecef[i] = _state(i, 0);
    // out << _var_state << endl << endl;
    // out.close();
    return true;
}

void CRtkekf::getDDApprox(sat* sats, MatrixXd state, double* b_pos, int* refsats, MatrixXd &w) {
    int nobs = sats->_nsats;
    int ref_prn = -1, sysflag = -1;
    int freqnum = _opt->_freqnum;
    int num = 0, ar_pos_ = 3;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = state(i, 0);
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1]._sat[isat]._isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((_opt->_navsys & SYS_ARRAY[i]) == sats[1]._sat[isat]._sys) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1]._sat[isat]._prn == ref_prn && sats[1]._sat[isat]._sys == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        double dist_base_ref = dist(ref_sat_b._pos, b_pos);
        double dist_base_com = dist(sats[0]._sat[isat]._pos, b_pos);
        double dist_rover_ref = dist(ref_sat_r._pos, rover_pos);
        double dist_rover_com = dist(sats[1]._sat[isat]._pos, rover_pos);
        double dist_cal = dist_rover_com - dist_rover_ref - dist_base_com + dist_base_ref;
        int sys = 0;
        for (int i_freq = 0; i_freq < MAXFREQ; ++i_freq) {
            if ((_opt->_freqtype & FREQ_ARRAY[i_freq]) == FREQ_ARRAY[i_freq]) {
                double freq = CPntbase::GetFreq(sysflag, FREQ_ARRAY[i_freq]);
                w(num ++, 0) = dist_cal;
                // states are: pos, GPS_Nf1, GPS_Nf2, BDS_NF1, BDS_NF2
                for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                    int total = 0;
                    for (int i_t = i_sys; i_t >= 0; -- i_t)
                        total += (_sysobs[i_t] - 1) * 2 * freqnum;
                    if (SYS_ARRAY[i_sys] == sysflag) {
                        w(num ++, 0) = dist_cal + state(ar_pos_ + i_freq * (_sysobs[i_sys] - 1), 0) * (VEL_LIGHT / freq);
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
    memset(_LeapSatsPos, 0, sizeof(int) * MAXOBS);
    _status = RTK_INIT;
    _opt = nullptr;
}

CRtkekf::CRtkekf(prcopt* opt) {
    if (!_opt) {
        _opt = new prcopt;
        memcpy(_opt, opt, sizeof(prcopt));
        _opt->_nbase = opt->_nbase; _opt->_nrover = opt->_nrover;
    }
    memset(_LeapSatsPos, 0, sizeof(int) * MAXOBS);
    _status = RTK_INIT;
}

MatrixXd CRtkekf::obssats(sat* sats_epoch, int nobs, int* refsats) {
    int nsats = sats_epoch->_nsats;
    MatrixXd usedsats(nobs, 3); usedsats.Zero();
    for (int i_sat = 0, i_row = 0; i_sat < nsats; ++i_sat) {
        if (!sats_epoch->_sat[i_sat]._isused) continue;
        usedsats(i_row, 0) = sats_epoch->_sat[i_sat]._sys;
        usedsats(i_row, 1) = sats_epoch->_sat[i_sat]._prn;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) 
            if (SYS_ARRAY[i_sys] == sats_epoch->_sat[i_sat]._sys && 
                refsats[i_sys] == sats_epoch->_sat[i_sat]._prn)
                usedsats(i_row, 2) = 1;
        i_row ++;
    }
    return usedsats;
}

void CRtkekf::updateStatus(MatrixXd obs, sat* sats_epoch, int nobs, int* refsats, res_t res, int* sysobs) {
    if (_status == RTK_INIT) {
        memcpy(_refsats, refsats, sizeof(int) * MAXSYS);
        memcpy(_sysobs, sysobs, sizeof(int) * MAXSYS);
        _state.resize((nobs - _opt->_nsys) * _opt->_freqnum + 3, 1);
        initvarsys(_state.row(), 1E-8);
        initvarstate(_state.row(), 15 * 15);
        initstate(sats_epoch, res._bpos_ecef, refsats, res);
        _state_trans.resize(_state.row(), _state.row()); 
        _state_trans.Identity();
        _nobs = nobs;
        _obssats = obssats(sats_epoch, nobs, refsats);
        _status = RTK_NORMAL;
    } else {
        MatrixXd obssat = obssats(sats_epoch, nobs, _refsats);
        if (checkObsChange(obssat, nullptr, nullptr))
            _status |= RTK_OBSCHANGE;
        for (int i = 0; i < MAXSYS; ++ i) 
            if (_refsats[i] != refsats[i]) _status |= RTK_REFCHANGE;
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
        _state(i, 0) = res._rpos_ecef[i];
    int ref_prn = -1, sysflag = -1, sys_pos = -1;
    int nobs = sats->_nsats;
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1]._sat[isat]._isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((_opt->_navsys & SYS_ARRAY[i]) == sats[1]._sat[isat]._sys) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i]; sys_pos = i;
                break;
            }
        }
        if(sats[1]._sat[isat]._prn == ref_prn && sats[1]._sat[isat]._sys == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        int nddobs = _sysobs[sys_pos] - 1;
        initambiguity(ref_sat_r, ref_sat_b, sats[1]._sat[isat], sats[0]._sat[isat], nddobs, ar_pos);
    }
}

void CRtkekf::initambiguity(sat_s ref_sat_r, sat_s ref_sat_b, sat_s sat_r, sat_s sat_b, int nddobs, int &ar_pos) {
    int num = 0;
    for (int i = 0; i < MAXFREQ; ++i) {
        if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
            double phase_obs = sat_r._obs->_L[i] - ref_sat_r._obs->_L[i]-
                sat_b._obs->_L[i] + ref_sat_b._obs->_L[i];
            double pseu_obs = sat_r._obs->_P[i] - ref_sat_r._obs->_P[i]- 
                sat_b._obs->_P[i] + ref_sat_b._obs->_P[i];
            double freq = CPntbase::GetFreq(sat_r._sys, FREQ_ARRAY[i]);
            double lambda = VEL_LIGHT / freq;
            if (num == 0)
                _state(ar_pos, 0) = -(pseu_obs - phase_obs * lambda) / lambda;
            else 
                _state(ar_pos + nddobs, 0) = -(pseu_obs - phase_obs * lambda) / lambda;
            ++ num;
        }
    }
    ar_pos += 1;
}

void CRtkekf::initvarsys(int dim, double var) {
    _var_sys.resize(dim, dim);
    for (int i = 0; i < dim; ++i) 
        for (int j = 0; j < dim; ++ j)
            if (i == j) _var_sys(i, j) = var;
            else        _var_sys(i, j) = 0;
}

void CRtkekf::initvarstate(int dim, double var) {
    _var_state.resize(dim, dim);
    for (int i = 0; i < dim; ++i) 
        for (int j = 0; j < dim; ++ j)
            if (i == j && i < 3) _var_state(i, j) = var;
            else if (i == j && i >= 3) _var_state(i, j) = var / 0.04;
            else        _var_state(i, j) = 0;
}

void CRtkekf::setopt(prcopt* opt) {
    if (!opt) return;
    _opt = opt;
}

CRtkekf::~CRtkekf() {
    _opt = nullptr;
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
        for (int j = 0; j < _obssats.row(); ++ j) {
            if (obssats_current(i, 0) == _obssats(j, 0) && 
                obssats_current(i, 1) == _obssats(j, 1)) {
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
    for (int i = 0; i < _obssats.row(); ++ i) {
        bool isfind = false;
        for (int j = 0; j < obssats_current.row(); ++ j) {
            if (obssats_current(j, 0) == _obssats(i, 0) && 
                obssats_current(j, 1) == _obssats(i, 1)) {
                isfind = true; break;
            }
        }
        if (isfind) continue;
        isdecrease = true;
        if (!decrease) continue; 
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys)
            if (SYS_ARRAY[i_sys] == int(_obssats(i, 0)))
                (*decrease)(sys_change[i_sys] ++, i_sys) = _obssats(i, 1);
    }
    return isdecrease;
}


void CRtkekf::statusFix(sat* sats_epoch, int nobs, int* refsats, int* obssys) {
    for (int i = 0; i < RTK_STATUS_NUM; ++i) {
        if ((_status & STATUS_ARRAY[i]) == STATUS_ARRAY[i]) 
            switch (STATUS_ARRAY[i]) {
            case RTK_LEAP:
                break;
            case RTK_REFCHANGE:
                refchange(refsats, obssys); _status &= 0x3B;
                memcpy(_refsats, refsats, sizeof(int) * MAXSYS);

                break;
            case RTK_OBSCHANGE: {
                MatrixXd obssats_epoch = obssats(sats_epoch, nobs, _refsats);
                MatrixXd increase;
                MatrixXd decrease;
                if (checkObsIncrease(obssats_epoch, &increase)) {
                    fixNumIncrease(sats_epoch, obssys, nobs, increase);
                    _nobs = nobs;
                }
                if (checkObsDecrease(obssats_epoch, &decrease)) {
                    fixNumDecrease(sats_epoch, obssys, nobs, decrease);
                    _nobs = nobs;
                }
                _status = RTK_NORMAL;
                break;
            }
            default:
                break;
            }
    }
    memcpy(_sysobs, obssys, sizeof(int) * MAXSYS);
    _obssats = obssats(sats_epoch, nobs, _refsats);
    _status = RTK_NORMAL;
}

void CRtkekf::fixNumIncrease(sat* sats, int* obssys, int nobs, MatrixXd increase) {
    int pos = -1, num_sys_obs = 0;
    int freqnum = _opt->_freqnum;
    bool isAfterRef = false;

    for (int i = 0; i < sats->_nsats; ++ i) {
        if (!sats->_sat[i]._isused) continue;
        bool isref = false;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (SYS_ARRAY[i_sys] == sats[1]._sat[i]._sys &&
                sats[1]._sat[i]._prn == _refsats[i_sys])
                isref = true;
        }
        if (isref) {
            isAfterRef = true;
            continue;
        } 
        pos += 1;
        if (!havesat(sats->_sat[i], increase)) continue;

        // construct transfer matrix
        MatrixXd T(_var_state.row() + freqnum, _var_state.col());
        MatrixXd T_info(_obssats.row() + 1, _obssats.row());  // transfer matrix info
        T.Zero(); T_info.Zero();
        for (int j = 0; j < 3; ++ j) T(j, j) = 1;
        int row = 3, col = 3;
        int sysflag = 0;
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys)
            if (sats->_sat[i]._sys == SYS_ARRAY[i_sys])
                sysflag = i_sys;
        int total_num = 0;
        for (int tmp = 0; tmp < sysflag; ++ tmp) {
            total_num += _sysobs[tmp] - 1; 
        }
        int pos_add = pos - total_num;
        // if (isAfterRef) pos_add -= 1;
        for (int i_part = 0; i_part < _opt->_nsys * _opt->_freqnum;  ++ i_part) {
            int sys_dd_num = _sysobs[int(i_part / freqnum)]; // contain new sats and ref sats
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
        _sysobs[sysflag] += 1;
        _state = T * _state;
        double ddobs[MAXSYS] = {0};
        // ofstream out("./debug_increase.txt");
        // auto test = obssats(sats, nobs, _refsats);
        // out << test << endl << endl;
        // out << _obssats << endl << endl;
        // out << increase << endl << endl;
        initstate_s(sats, sats->_sat[i]._prn, sats[1]._sat[i]._sys, ddobs);
        // out << _state << endl << endl;
        // out << ddobs[0] << "  " << ddobs[1] << endl << endl;
        for (int j = 0, num = 0; j < _state.row(); ++ j)
            if (_state(j, 0) == 0)
                _state(j, 0) = ddobs[num ++];

        _var_state = T * _var_state * T.transpose();
        _state_trans.resize(_state.row(), _state.row()); _state_trans.Identity();
        _var_sys.resize(_state.row(), _state.row()); _var_sys.Zero();
        // out << T << endl << endl;
        // out << _state << endl << endl;
        // out << _var_state << endl << endl;
        // out.close();
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
        _obssats = T_info * _obssats;
        for (int i_row = 0; i_row < _obssats.row(); ++ i_row)
            if (_obssats(i_row, 0) == 0) {
                _obssats(i_row, 0) = sats->_sat[i]._sys;
                _obssats(i_row, 1) = sats->_sat[i]._prn;
            }
        for (int i_row = 0; i_row < _var_state.row(); ++i_row) {
            if (abs(_var_state(i_row, i_row)) <= 1e-8)
                _var_state(i_row, i_row) = INIT_STATE_VAR * INIT_STATE_VAR / 0.19 * 100;
        }
        pos = tmp;  // for any number of satellites increase
    }
    _nobs = nobs;
}

void CRtkekf::fixNumDecrease(sat* sats, int* obssys, int nobs, MatrixXd decrease) {
    int num_sys_obs = 0;
    int freqnum = _opt->_freqnum;
    bool isafterRef = false;
    int ref_pos[MAXSYS];
    for (int i = 0; i < MAXSYS; i ++) ref_pos[i] = -1;
    for (int i = 0; i < _obssats.row(); ++ i) {
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (_refsats[i_sys] == _obssats(i, 1) && 
                SYS_ARRAY[i_sys] == _obssats(i, 0) && 
                _obssats(i, 2) == 1)
                ref_pos[i_sys] = i;
        }
    }
    // ofstream out ("./debug_decrease.txt");
    // out << decrease << endl << endl;
    // out << _obssats << endl << endl;
    for (int i_row_c = 0; i_row_c < decrease.row(); ++ i_row_c) {
        for (int i_row_o = 0; i_row_o < _obssats.row(); ++ i_row_o) {
            for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
                if (SYS_ARRAY[i_sys] != _obssats(i_row_o, 0) ||
                    _obssats(i_row_o, 1) != decrease(i_row_c, i_sys))
                    continue;
                
                MatrixXd T_trans(_state.row() - freqnum, _state.row()); T_trans.Zero();
                for (int i = 0; i < 3; ++ i) T_trans(i, i) = 1;
                int pos = i_row_o, total_obs = 0;
                int row = 3, col = 3;
                for (int i = 0; i < i_sys; ++ i) total_obs += _sysobs[i] - 1;
                pos -= total_obs == 0 ? total_obs : (total_obs + 1);
                ref_pos[i_sys] -= (total_obs == 0 ? total_obs : total_obs + 1);
                if (pos >= ref_pos[i_sys]) {
                    pos -= 1; isafterRef = true;
                }
                for (int i_part = 0; i_part < _opt->_freqnum * _opt->_nsys; ++ i_part) {
                    if (i_part == i_sys * 2 || i_part - 1 == i_sys * 2) {
                        int ddobs = _sysobs[i_part / freqnum] - 1;
                        MatrixXd T_part (ddobs - 1, ddobs); T_part.Zero();
                        for (int i = 0; i < T_part.row(); ++ i) 
                            if (i >= pos) T_part(i, i + 1) = 1;
                            else T_part(i, i) = 1;
                        T_trans.block(row, col, T_part);
                        row += T_part.row(); col += T_part.col();
                    } else if (i_part != i_sys * 2 || i_part - 1 != i_sys * 2) {
                        int ddobs = _sysobs[i_part / freqnum] - 1;
                        MatrixXd T_part (ddobs, ddobs); T_part.Identity();
                        T_trans.block(row, col, T_part);
                        row += T_part.row(); col += T_part.col();
                    }
                }
                // out << T_trans << endl << endl;
                _sysobs[i_sys] -= 1;
                // out << _state << endl << endl;
                _state = T_trans * _state;
                _var_state = T_trans * _var_state * T_trans.transpose();
                // out << _state << endl << endl;
                if (isafterRef) pos += 1;
                pos += (total_obs == 0 ? total_obs : total_obs + 1);
                MatrixXd T_info(_obssats.row() - 1, _obssats.row()); T_info.Zero();
                for (int i = 0; i < T_info.row(); ++ i) {
                    if (i < pos) T_info(i, i) = 1;
                    else T_info(i, i + 1) = 1;
                }
                _obssats = T_info * _obssats;
                // out << _obssats << endl << endl;
                _state_trans.resize(_state.row(), _state.row()); _state_trans.Identity();
                _var_sys.resize(_state.row(), _state.row()); _var_sys.Zero();
                _nobs = nobs;
            }
        }
    }
}


bool CRtkekf::havesat(sat_s sat, MatrixXd changeObs) {
    for (int i_row = 0; i_row < MAX_OBSCHANGE; ++ i_row) {
        for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
            if (SYS_ARRAY[i_sys] == sat._sys && 
                sat._prn == changeObs(i_row, i_sys)) 
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
    int freqnum = _opt->_freqnum;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i]);

        for (int i_row = 0; i_row < _obssats.row(); ++ i_row) 
            if (_obssats(i_row, 0) == SYS_ARRAY[i_sys] && 
                _obssats(i_row, 1) == refsats[i_sys] && 
                _obssats(i_row, 2) != 1) {
                // position of L1-ambiguity 
                new_ref_position[i_sys] = i_row - total; 
            }
        
        for (int i_row = 0; i_row < _obssats.row(); ++ i_row) 
            if (_obssats(i_row, 0) == SYS_ARRAY[i_sys] && 
                _obssats(i_row, 2) == 1)
                if (i_row <= new_ref_position[i_sys] + total)
                    new_ref_position[i_sys] -= 1;
    }
    
    MatrixXd trans(_state.row(), _state.row()); trans.Zero();
    for (int i = 0; i < 3; ++ i) trans(i, i) = 1;

    int row = 3;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += (obssys[i] - 1) * freqnum;
        int ddobs = obssys[i_sys] - 1;
        if (_refsats[i_sys] == refsats[i_sys] || new_ref_position[i_sys] == -1) {
            MatrixXd identity(ddobs, ddobs);
            identity.Identity();
            for (int i_freq = 0; i_freq < freqnum; ++ i_freq)
                trans.block(row + total + i_freq * (obssys[i_sys] - 1), row + total + i_freq * (obssys[i_sys] - 1), identity);
            continue;
        }
        MatrixXd trans_part(ddobs, ddobs);

        for (int ifreq = 0; ifreq < _opt->_freqnum; ++ ifreq) {
            trans_part.Zero();
            for (int i = 0; i < trans_part.row(); ++ i) {
                trans_part(i, new_ref_position[i_sys]) = -1;
                if (i == new_ref_position[i_sys]) continue;
                trans_part(i, i) = 1;
            }
            trans.block(row + total + ddobs * ifreq, row + total + ddobs * ifreq, trans_part);
        }
        
    }
    _state = trans * _state;
    _var_state = trans * _var_state * trans.transpose();
    sortRefChange(new_ref_position, refsats, obssys);
}

void CRtkekf::sortRefChange(int* new_ref_position, int* refsats, int* obssys) {
    int old_ref_position[MAXSYS];
    int freqnum = _opt->_freqnum;
    for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) {
        if (new_ref_position[i_sys] == -1) continue;
        int total = 0;
        for (int i = 0; i < i_sys; i ++) total += obssys[i];

        for (int i_row = 0; i_row < _obssats.row(); ++ i_row) 
            if (_obssats(i_row, 0) == SYS_ARRAY[i_sys] && 
                _obssats(i_row, 1) == _refsats[i_sys]) {
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
            // for the part will be sorted to _obssats position
            MatrixXd part_state(ddobs, _state.col());
            MatrixXd part_var_state(ddobs, ddobs);

            for (int i = 0; i < ddobs; ++ i)
                part_state(i, 0) = _state(row + i, 0);
            for (int i = 0; i < ddobs; ++ i)
                for (int j = 0; j < part_var_state.col(); ++ j)
                    part_var_state(i, j) = _var_state(row + i, row + j);
            
            MatrixXd temp_state(1, _state.col());
            int new_ref_pos = new_ref_position[i_sys];
            int old_ref_pos = old_ref_position[i_sys];
            temp_state(0, 0) = part_state(new_ref_pos, 0);
            MatrixXd T_part(part_state.row(), part_state.row()); T_part.Zero();
            MatrixXd T(_var_state.row(), _var_state.col());

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
            _var_state = T * _var_state * T.transpose();
            _state.block(row, 0, part_state);
        }
    }

}

void CRtkekf::initstate_s(sat* sats, int prn, int sys, double* ddobs) {
    int nobs = sats->_nsats;
    int ref_prn = -1, sysflag = -1;
    int freqnum = _opt->_freqnum;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = _state(i, 0);
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1]._sat[isat]._isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((_opt->_navsys & SYS_ARRAY[i]) == sats[1]._sat[isat]._sys) {
                ref_prn = _refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if (sats->_sat[isat]._sys != sys || sats->_sat[isat]._prn != prn)
            continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        sat_s sat_r = sats[1]._sat[isat], sat_b = sats[0]._sat[isat];
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sat_r._obs->_L[i] - ref_sat_r._obs->_L[i]-
                    sat_b._obs->_L[i] + ref_sat_b._obs->_L[i];
                double pseu_obs = sat_r._obs->_P[i] - ref_sat_r._obs->_P[i]- 
                    sat_b._obs->_P[i] + ref_sat_b._obs->_P[i];
                double freq = CPntbase::GetFreq(sat_r._sys, FREQ_ARRAY[i]);
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
    if (_opt->_freqnum < 2) return; // at least dual-freq
    for (int i_site = 0; i_site < _opt->_sitenum; ++ i_site) {
        int pos = -1;
        bool isafterref = false;
        for (int i_sat = 0; i_sat < sats_epoch[i_site]._nsats; ++ i_sat) {
            if (!sats_epoch[i_site]._sat[i_sat]._isused) continue;
            pos += 1;

            // gf combination
            double gf = 0;
            for (int i_freq = 0; i_freq < MAXFREQ; ++ i_freq) {
                if ((_opt->_freqtype & FREQ_ARRAY[i_freq]) == FREQ_ARRAY[i_freq]) {
                    double freq = CPntbase::GetFreq(sats_epoch[i_site]._sat[i_sat]._sys, FREQ_ARRAY[i_freq]);
                    double lambda = VEL_LIGHT / freq;
                    gf += sats_epoch[i_site]._sat[i_sat]._obs->_L[i_freq] * lambda * pow(-1, i_freq);
                }
            }
            sats_epoch[i_site]._sat[i_sat]._gf = gf;
            if (!sat_last_epoch) continue; // not initialized 
            int sysflag = 0;
            for (int i_sys = 0; i_sys < MAXSYS; ++ i_sys) 
                if (SYS_ARRAY[i_sys] == sats_epoch[i_site]._sat[i_sat]._sys)
                    sysflag = i_sys;
            if (sats_epoch[i_site]._sat[i_sat]._prn == _refsats[sysflag]) isafterref = true;
            int state_pos = pos;
            int ddobs = 0, freqnum = _opt->_freqnum;
            for (int i_sys = 0; i_sys < sysflag; ++ i_sys) ddobs += (_sysobs[i_sys] - 1);
            state_pos -= ddobs == 0 ? 0 : (ddobs + 1); // minus ref sat
            state_pos += freqnum * ddobs;   // restore the true position in state
            if (isafterref) state_pos -= 1;
            state_pos += 3;
            double last_gf = searchgf(sat_last_epoch[i_site], sats_epoch[i_site]._sat[i_sat]);
            if (last_gf == -1) continue;

            if (abs(gf - last_gf) >= THRESHOLD) { // cycle slip occurred, reset state
                // ofstream out("./debug_slip.txt");
                // out << _obssats << endl << endl;
                // out.close();
                auto ref_sat_b = findRef(sats_epoch[0], sats_epoch[0]._sat[i_sat]._sys, _refsats[sysflag]);
                auto ref_sat_r = findRef(sats_epoch[1], sats_epoch[1]._sat[i_sat]._sys, _refsats[sysflag]);
                auto sat_b = sats_epoch[0]._sat[i_sat];
                auto sat_r = sats_epoch[1]._sat[i_sat];
                initambiguity(ref_sat_r, ref_sat_b, sat_r, sat_b, _sysobs[sysflag] - 1, state_pos);
                state_pos -= 1;
                for (int i = 0; i < _var_state.row(); ++ i) 
                    if (i != state_pos) _var_state(state_pos, i) = 0;
                    else _var_state(state_pos, state_pos) = 30 * 30;
                state_pos += _sysobs[sysflag] - 1;

                for (int i = 0; i < _var_state.row(); ++ i) 
                    if (i != state_pos) _var_state(state_pos, i) = 0;
                    else _var_state(state_pos, state_pos) = 30 * 30;
            }
        } // each satellite loop
    } // each site loop

    if (!sat_last_epoch) sat_last_epoch = new sat[_opt->_sitenum];
    memcpy(sat_last_epoch, sats_epoch, sizeof(sat) * _opt->_sitenum);
}

double CRtkekf::searchgf(sat sat_last, sat_s obj) {
    for (int isat = 0; isat < sat_last._nsats; ++ isat)
        if (sat_last._sat[isat]._isused && sat_last._sat[isat]._sys == obj._sys
            && sat_last._sat[isat]._prn == obj._prn)
            return sat_last._sat[isat]._gf;
    return -1;
}

int CRtkekf::obsnumber(sat* sats) {
    int nsites = _opt->_sitenum;
    int nobs = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < sats->_nsats; ++ i_sat) {
        if (!sats[0]._sat[i_sat]._isused || !sats[1]._sat[i_sat]._isused)
            sats[0]._sat[i_sat]._isused = sats[1]._sat[i_sat]._isused = false;
        else nobs += 1;
    }
    return nobs;
}

void CRtkekf::getDesignDim(sat sats, int nobs, int &row, int &col) {
    int nsys = _opt->_nsys, nfreq = _opt->_freqnum;
    row = col = 0;    
    row = (nobs - _opt->_nsys) * nfreq * 2;
    col = (nobs - _opt->_nsys) * nfreq + 3;
}

void CRtkekf::getddobs(sat* sats, double* sitepos, int* refsats, int* sysobs, res_t res, MatrixXd &w) {
    int nobs = sats->_nsats;
    int ref_prn = -1, sysflag = -1;
    int num = 0;
    double rover_pos[3] = {0};
    for (int i = 0; i < 3; ++i) rover_pos[i] = res._rpos_ecef[i];
    for (int isat = 0; isat < nobs; ++isat) {
        if(!sats[1]._sat[isat]._isused)  continue;
        for (int i = 0; i < MAXSYS; ++i) {
            if ((_opt->_navsys & SYS_ARRAY[i]) == sats[1]._sat[isat]._sys) {
                ref_prn = refsats[i]; sysflag = SYS_ARRAY[i];
                break;
            }
        }
        if(sats[1]._sat[isat]._prn == ref_prn && sats[1]._sat[isat]._sys == sysflag) continue;
        sat_s ref_sat_r = findRef(sats[1], sysflag, ref_prn);
        sat_s ref_sat_b = findRef(sats[0], sysflag, ref_prn);
        for (int i = 0; i < MAXFREQ; ++i) {
            if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i]) {
                double phase_obs = sats[1]._sat[isat]._obs->_L[i] - ref_sat_r._obs->_L[i]-
                    sats[0]._sat[isat]._obs->_L[i] + ref_sat_b._obs->_L[i];
                double pseu_obs = sats[1]._sat[isat]._obs->_P[i] - ref_sat_r._obs->_P[i]- 
                    sats[0]._sat[isat]._obs->_P[i] + ref_sat_b._obs->_P[i];
                double freq = CPntbase::GetFreq(sysflag, FREQ_ARRAY[i]);
                w(num ++, 0) = pseu_obs;
                w(num ++, 0) = phase_obs * (VEL_LIGHT / freq);
            }
        }
    }
}

MatrixXd CRtkekf::GetState() {
    return _state;
}

MatrixXd CRtkekf::GetVar() {
    return _var_state;
}

double CRtkekf::GetInternalSigma() {
    return 1.0;
}
