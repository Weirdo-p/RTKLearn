#include "navigation/pnt/pntrtk.h"
#include "navigation/timemodule.h"
#include "navigation/coors.h"
#include <string.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <chrono>
#include "navigation/optimal/rtkekf.h"
#include "navigation/ambiguity/lambda.h"

using namespace chrono;

CPntrtk::CPntrtk() {
    _opt = nullptr; _optimizer = nullptr;
}

CPntrtk::CPntrtk(prcopt opt) {
    if (_opt) {
        memcpy(_opt, &opt, sizeof(prcopt));
    } else {
        _opt = new prcopt;
        memcpy(_opt, &opt, sizeof(prcopt));
    }
    if (!_optimizer) {
        if (_opt->_proctype == PROC_KF)
            _optimizer = new CRtkekf(_opt);
        if (_opt->_proctype == PROC_LS)
            _optimizer = new CLeastsq(_opt);
    }
    _spprunner = new CPntspp(*_opt);
}

CPntrtk::~CPntrtk() {
    if (_spprunner) 
        delete _spprunner;
    _spprunner = nullptr;
}

int CPntrtk::process() {
    if (_opt->_sitenum < 2) return -1;
    sat* sats_epoch = new sat[_opt->_sitenum];
    if (_opt->_base[0] != 0) memcpy(_res->_bpos_ecef, _opt->_base, sizeof(double) * 3);
    if (_opt->_rover[0] != 0) memcpy(_res->_rpos_ecef, _opt->_rover, sizeof(double) * 3);
    
    while(inputobs(sats_epoch)) {
        auto start = system_clock::now();
        satclk(sats_epoch);
        satpos(sats_epoch);
        earthRotateFix(sats_epoch);
        satvel(sats_epoch);
        relativeeffect(sats_epoch);
        _spprunner->spp(sats_epoch);
        auto res = _spprunner->getRes();
        memcpy(_res->_rpos_ecef, res->_rpos_ecef, sizeof(double) * 3);
        satazel(_res->_rpos_ecef, sats_epoch[0]);
        rtk(sats_epoch);
        auto end = system_clock::now();
        auto cost = duration_cast<microseconds> (end - start);
        cout << setw(4) << sats_epoch->_sat->_obs->_time._Week << " " << setw(9) << fixed << setprecision(6) << sats_epoch->_sat->_obs->_time._Sow;
        cout<< " TIME COST: " << double (cost.count()) * microseconds::period::num / microseconds::period::den << " seconds" << endl;
        memset(sats_epoch, 0, sizeof(sat) * _opt->_sitenum);
    }
}

void CPntrtk::outsol(Sattime time, int nobs) {
    Commontime com;
    unsigned short doy = 0;
    GPST2Common(time, com, 0); Common2Doy(com, doy);
    string file_name = "./results/result_" + to_string(doy) + ".txt";
    ofstream out(file_name, ios::app);
    out << setfill(' ');
    out << setw(4) << time._Week << " " << setw(18) << fixed << setprecision(6) << time._Sow;
    out << " " << setw(18) << fixed << setprecision(6) << _res->_rpos_ecef[0] <<
            " " << setw(18) << fixed << setprecision(6) << _res->_rpos_ecef[1] <<
            " " << setw(18) << fixed << setprecision(6) << _res->_rpos_ecef[2] << 
            " " << setw(9) << fixed << setprecision(6) << _res->_enu[0] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_enu[1] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_enu[2] << 
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[0] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[1] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[2] <<
            " " << setw(2) << fixed << _res->_ambi_flag << 
            " " << setw(5) << fixed << setprecision(3) << _res->_ratio << 
            " " << setw(5) << fixed << setprecision(3) << _res->_rdop;
    out << setw(3) << fixed << nobs << endl;
    out.close();
}

int CPntrtk::inputobs(sat* sats) {
    // record observations' postion
    static int base_pos = 0, rover_pos = 0;
    int nsites = _opt->_sitenum;
    obs* obss = _rnx->GetObs();
    nav* navs = _rnx->GetEph();
    obs* base = &obss[0], *rover = &obss[1];
    double diff = Sattimediff(base->_obs[base_pos]._time, rover->_obs[rover_pos]._time);
    // searching observations
    while(abs(diff) >= 1e-3) {
        if (base_pos >= base->_obsnum || rover_pos >= rover->_obsnum)
            break;
        if (diff < 0)   base_pos += 1;
        else            rover_pos += 1;
        diff = Sattimediff(base->_obs[base_pos]._time, rover->_obs[rover_pos]._time);
    }
    // searching ephemeris
    if (base_pos >= base->_obsnum || rover_pos >= rover->_obsnum)
        return 0;
    setobs(obss, base_pos, rover_pos, sats);
    setnav(navs, sats);
    obss = nullptr; navs = nullptr; base = nullptr; rover = nullptr;
    return 1;
}

int CPntrtk::setobs(obs* obss, int &base_pos, int &rover_pos, sat* sats) {
    obs* base = &obss[0], *rover = &obss[1];
    int nsats = 0;
    Sattime time = base->_obs[base_pos]._time;
    int i_base = base_pos;
    #pragma omp parallel for
    for (; i_base < base->_obsnum; ++ i_base) {
        if (base->_obs[i_base]._time != time) break;
        for (int i_rover = rover_pos; i_rover < rover->_obsnum; ++ i_rover) {
            if (rover->_obs[i_rover]._time != time) break;
            if (base->_obs[i_base]._sys != rover->_obs[i_rover]._sys || 
                base->_obs[i_base]._sat != rover->_obs[i_rover]._sat)
                continue;
            sats[0]._sat[nsats]._prn = sats[1]._sat[nsats]._prn = rover->_obs[i_rover]._sat;
            sats[0]._sat[nsats]._sys = sats[1]._sat[nsats]._sys = rover->_obs[i_rover]._sys;
            sats[0]._sat[nsats]._obs = &base->_obs[i_base];
            sats[1]._sat[nsats]._obs = &rover->_obs[i_rover];
            nsats += 1;
        }
    }
    sats[0]._nsats = sats[1]._nsats = nsats;
    base_pos = i_base;
    base = nullptr; rover = nullptr;
    return 1;
}

int CPntrtk::setnav(nav* navs, sat* sats) {
    for (int i = 0; i < 2; ++i) {
        #pragma omp parallel for
        for (int isat = 0; isat < sats[i]._nsats; ++isat) {
            int sys = sats[i]._sat[isat]._sys;
            int prn = sats[i]._sat[isat]._prn;
            Sattime time = sats[i]._sat[isat]._obs->_time;
            sats[i]._sat[isat]._eph = searchnav(time, sys, prn, navs);
        }
    }
    return 1;
}

bool CPntrtk::rtk(sat* sats_epoch) {
    int MAXITER = 10, iter = 0, issuccess = 0;
    int REFSATS[MAXSYS] = { -1 };
    int nobs = 0;
    for (int i = 0; i < _opt->_sitenum; ++i)
        excludesats(sats_epoch[i]);  // observation count by rover
    
    nobs = obsnumber(sats_epoch);
    if (sats_epoch->_sat->_obs->_time._Sow == 551406)
        cout << "test" << endl;
    _optimizer->optimize(sats_epoch, *_res);
    evaluate();
    fixambi();
    XYZ2NEU(_res->_bpos_ecef, _res->_rpos_ecef, WGS84, _res->_enu);
    outsol(sats_epoch->_sat->_obs->_time, nobs);
}

int CPntrtk::obsnumber(sat* sats) {
    int nobs = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < sats->_nsats; ++ i_sat) {
        if (!sats[0]._sat[i_sat]._isused || !sats[1]._sat[i_sat]._isused)
            sats[0]._sat[i_sat]._isused = sats[1]._sat[i_sat]._isused = false;
        else nobs += 1;
    }
    return nobs;
}

int CPntrtk::excludesats(sat &sat) {
    int satnum = sat._nsats;
    double cutoff = _opt->_elecutoff;
    for(int isat = 0; isat < satnum; ++ isat) {
        if(sat._sat[isat]._elev != 0 && sat._sat[isat]._elev < _opt->_elecutoff) {
            sat._sat[isat]._isused = false; continue;
        }
        if ((sat._sat[isat]._prn <= 5) && sat._sat[isat]._sys == SYS_BDS) {
            sat._sat[isat]._isused = false; continue;
        }
        if ((sat._sat[isat]._prn >= 59) && sat._sat[isat]._sys == SYS_BDS) {
            sat._sat[isat]._isused = false; continue;
        }

        else {
            bool isobs = true;
            for (int i = 0; i < MAXFREQ; ++i)
                if ((_opt->_freqtype & FREQ_ARRAY[i]) == FREQ_ARRAY[i] && 
                    abs(sat._sat[isat]._obs->_P[i]) <= 1e-6) {
                        isobs = false; break;
                    }
            if (isobs) {
                sat._sat[isat]._isused = true;
            }
            else sat._sat[isat]._isused = false;
        }
        if (!sat._sat[isat]._obs || !sat._sat[isat]._eph)
            sat._sat[isat]._isused = false;
    }
}

void CPntrtk::evaluate() {
    MatrixXd var_state = _optimizer->GetVar();
    MatrixXd Q = var_state.block<0, 0>(3, 3);
    double sigma = _optimizer->GetInternalSigma();

    double station_blh[3] = {0};
    XYZ2BLH(_res->_bpos_ecef, WGS84, station_blh);
    double B = station_blh[0], L = station_blh[1];
    MatrixXd rotation(3, 3);
    rotation(0, 0) = -sin(B) * cos(L); rotation(0, 1) = -sin(B) * sin(L); rotation(0, 2) = cos(B);
    rotation(1, 0) = -sin(L); rotation(1, 1) = cos(L); rotation(1, 2) = 0;
    rotation(2, 0) = cos(B) * cos(L); rotation(2, 1) = cos(B) * sin(L); rotation(2, 2) = sin(B);

    // project to neu coordinate
    MatrixXd Q_neu = rotation * Q * rotation.transpose();   
    Q_neu = Q_neu * sigma;
    for (int i = 0; i < 3; ++ i)
        _res->_sigma_neu[i] = sqrt(Q_neu(i, i));
    
    // RDOP
    double rdop = 0;
    for (int i = 0; i < var_state.row(); ++ i) rdop += var_state(i, i) * sigma;
    _res->_rdop = sqrt(rdop);
}
