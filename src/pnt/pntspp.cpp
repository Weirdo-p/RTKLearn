#include "navigation/pnt/pntspp.h"
#include "navigation/atmosphere.h"
#include "navigation/coors.h"
#include "navigation/timemodule.h"
#include <iomanip>
#include <fstream>
#include <chrono>
using namespace chrono;

CPntspp::CPntspp() {

}

CPntspp::CPntspp(prcopt opt) {
    if (_opt)
        memcpy(_opt, &opt, sizeof(opt));
    else {
        _opt = new prcopt;
        memcpy(_opt, &opt, sizeof(opt));
    }
    _optimizer = new CLeastsq(&opt);
}

CPntspp::~CPntspp() {
    if(_optimizer) delete _optimizer;
    _optimizer = nullptr;
}


int CPntspp::process() {
    if (_opt->_sitenum < 1) return -1;
    sat* sats_epoch = new sat;
    if (_opt->_rover[0] != 0) memcpy(_res->_rpos_ecef, _opt->_rover, sizeof(double) * 3);
    
    while(inputobs(sats_epoch)) {
        auto start = system_clock::now();
        satclk(sats_epoch);
        satpos(sats_epoch);
        earthRotateFix(sats_epoch);
        satvel(sats_epoch);
        relativeeffect(sats_epoch);
        bool issuccess = true;
        for (int i_site = 0; i_site < _opt->_sitenum; ++ i_site) {
            if (spp_site(i_site, &sats_epoch[i_site]) == -1) issuccess = false;
            avd(i_site, &sats_epoch[i_site]);
        }
        if (!issuccess) continue;
        int nobs = excludesats(*sats_epoch);
        outsol(sats_epoch->_sat->_obs->_time, nobs);
        auto res = getRes();
        memcpy(_res->_rpos_ecef, res->_rpos_ecef, sizeof(double) * 3);
        auto end = system_clock::now();
        auto cost = duration_cast<microseconds> (end - start);
        cout << setw(4) << sats_epoch->_sat->_obs->_time._Week << " " << setw(9) << fixed << setprecision(6) << sats_epoch->_sat->_obs->_time._Sow;
        cout<< " TIME COST: " << double (cost.count()) * microseconds::period::num / microseconds::period::den << " seconds" << endl;
        memset(sats_epoch, 0, sizeof(sat) * _opt->_sitenum);
    }
}


int CPntspp::inputobs(sat* sats) {
    static int rover_pos = 0;
    int nsites = _opt->_sitenum;
    obs* obss = _rnx->GetObs();
    nav* navs = _rnx->GetEph();
    obs* rover = &obss[0];

    setobs(obss, rover_pos, sats);
    setnav(navs, sats);
    if (rover_pos >= rover->_obsnum)
        return 0;
    obss = nullptr; navs = nullptr; rover = nullptr;
    return 1;
}

int CPntspp::setobs(obs* obss, int &rover_pos, sat* sats) {
    obs* rover = &obss[0];
    int nsats = 0;
    Sattime time = rover->_obs[rover_pos]._time;
    int i_rover = rover_pos;
    #pragma omp parallel for
    for (; i_rover < rover->_obsnum; ++ i_rover) {
        if (rover->_obs[i_rover]._time != time) break;
        sats[0]._sat[nsats]._prn = sats[1]._sat[nsats]._prn = rover->_obs[i_rover]._sat;
        sats[0]._sat[nsats]._sys = sats[1]._sat[nsats]._sys = rover->_obs[i_rover]._sys;
        sats[0]._sat[nsats]._obs = &rover->_obs[i_rover];
        nsats += 1;
    }

    sats[0]._nsats = sats[1]._nsats = nsats;
    rover_pos = i_rover; rover = nullptr;
    return 1;
}

int CPntspp::setnav(nav* navs, sat* sats) {
    for (int isat = 0; isat < sats[0]._nsats; ++isat) {
        int sys = sats[0]._sat[isat]._sys;
        int prn = sats[0]._sat[isat]._prn;
        Sattime time = sats[0]._sat[isat]._obs->_time;
        sats[0]._sat[isat]._eph = searchnav(time, sys, prn, navs);
    }
    return 1;
}

int CPntspp::spp(sat* sats) {
    int sitenum = _opt->_sitenum;
    int start = 0;
    if (_opt->_mode == MODE_SINGLE) start = 0;
    if (_opt->_base[0] != 0 && _opt->_mode == MODE_RTK) start = 1;
    for (int i_site = start; i_site < _opt->_sitenum; ++ i_site) {
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
            H = _res->_bpos_blh[2];
        } else {
            H = _res->_rpos_blh[2];
        }
        int nsys = usesys(*sats);
        int cols = 4;
        for (int i = 0; i < nsys - 1; ++i)
            cols += 1;
        if (nobs < cols) {
            cout << "FATAL: no sufficient observations" << endl;
            return -1;
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
        if(!_optimizer->optimize(B, P, w, nobs)) {
            // memset(sitepos_ecef, 0, sizeof(double) * 3);
            // cout << w << endl;
            // cout << B << endl;
            // return -1;
        }
        
        MatrixXd x = _optimizer->GetState();
        auto site_pos = x.block<0, 0>(3, 1);
        for (int i = 0; i < x.row(); ++i) pos(i, 0) += x(i, 0);
        setres(sitepos_ecef, sitepos_blh, type, pos);
        if (site_pos.norm() < 1E-6)
            break;
    }
    memcpy(_res->_rpos_ecef, sitepos_ecef, sizeof(double) * 3);
    XYZ2BLH(sitepos_ecef, type, _res->_rpos_blh);
    evaluate(); // TODO:
    sitepos_ecef = nullptr; sitepos_blh = nullptr;
    return 0;
}

void CPntspp::GetDesign(sat sats, double* sitepos, MatrixXd &B) {
    int satnum = sats._nsats, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats._sat[i_sat]._isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats._sat[i_sat]._pos[i];
            dist += coeff[i] * coeff[i];
        }
        int p_pos = 0;  // used pseudorange position
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats._sat[i_sat]._sys & SYS_ARRAY[i]) == SYS_ARRAY[i])
                p_pos = i;
        if (_opt->_nsys == 1) p_pos = 0;
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
    int satnum = sats._nsats, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats._sat[i_sat]._isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats._sat[i_sat]._pos[i];
            dist += coeff[i] * coeff[i];
        }
        int p_pos = 0;
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats._sat[i_sat]._sys & SYS_ARRAY[i]) == SYS_ARRAY[i])
                p_pos = i;
        if (_opt->_nsys == 1) p_pos = 0;
        dist = sqrt(dist);
        w(num++, 0) = -dist + sats._sat[i_sat]._clk[0] * VEL_LIGHT - pos(3 + p_pos, 0);
    }
    if (_opt->_freqnum != 1) GetLC(sats, w);
    else GetNonCombine(sats, w);
    double BLH[3];
    if (*sitepos == 0) return;
    XYZ2BLH(sitepos, Ellipsoid(WGS84), BLH);
    if (abs(BLH[2]) > 1E4) 
        return;
    tropfix(sats, w, BLH[2]);
}

res_t* CPntspp::getRes() {
    return _res;
}

int CPntspp::excludesats(sat &sat) {
    int satnum = sat._nsats;
    int n = 0;
    double cutoff = _opt->_elecutoff;
    for(int isat = 0; isat < satnum; ++ isat) {
        if(sat._sat[isat]._elev != 0 && sat._sat[isat]._elev < _opt->_elecutoff) {
            sat._sat[isat]._isused = false; continue;
        }
        if (!sat._sat[isat]._eph) {
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
                n += 1;
                sat._sat[isat]._isused = true;
            }
            else sat._sat[isat]._isused = false;
        }
    }

    return n;
}

void CPntspp::outsol(Sattime time, int nobs) {
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
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[0] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[1] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_neu[2] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_vel[0] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_vel[1] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_vel[2] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_vel[0] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_vel[1] <<
            " " << setw(9) << fixed << setprecision(6) << _res->_sigma_vel[2] <<
            " " << setw(2) << fixed << _res->_ambi_flag << 
            " " << setw(5) << fixed << setprecision(3) << _res->_rdop;
    out << setw(3) << fixed << nobs << endl;
    out.close();
}

void CPntspp::evaluate() {
    MatrixXd var_state = _optimizer->GetVar();
    MatrixXd Q = var_state.block<0, 0>(3, 3);
    double sigma = _optimizer->GetInternalSigma();

    double station_blh[3] = {0};
    XYZ2BLH(_res->_rpos_ecef, WGS84, station_blh);
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
    for (int i = 0; i < 3; ++ i) rdop += Q(i, i);
    _res->_rdop = sqrt(rdop);
    _res->_ambi_flag = SPP_SOLU;
}

void CPntspp::avd(int i_site, sat* sats) {
    int MAXITER = 10;
    double* sitepos_ecef, *sitepos_blh;

    selectPos(&sitepos_ecef, &sitepos_blh, i_site);
    int nobs = excludesats(*sats);
    MatrixXd pos;
    Ellipsoid type(WGS84);
    int nsys = usesys(*sats);
    int cols = 4;
    for (int i = 0; i < nsys - 1; ++i)
        cols += 1;
    if (nobs < cols) {
        cout << "FATAL: no sufficient observations" << endl;
        return;
    }
    pos.resize(cols, 1); pos.Zero();
    MatrixXd B(nobs, cols), P(nobs, nobs);
    P.Identity(); B.Zero();
    MatrixXd w(nobs, 1);
    GetAvdl(*sats, sitepos_ecef, pos, w);
    // cout << w << endl << endl;
    GetAvdDesign(*sats, sitepos_ecef, B);
    // cout << B << endl << endl;
    // cout << pos << endl << endl;
    if(!_optimizer->optimizeAVD(B, P, w, nobs)) {
        // pos.Zero();
        // cout << w << endl;
        // cout << B << endl;
    }
    MatrixXd x = _optimizer->GetVel();
    for (int i = 0; i < x.row(); ++i) _res->_vel[i] = x(i, 0);

    evaluateAVD();
    // if(_opt->_mode == MODE_SINGLE)
    //     outsol(sats->_sat->_obs->_time, nobs);
    sitepos_ecef = nullptr; sitepos_blh = nullptr;
}

int CPntspp::GetAvdDesign(sat sats, double* sitepos, MatrixXd &B) {
    int satnum = sats._nsats, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats._sat[i_sat]._isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats._sat[i_sat]._pos[i];
            dist += coeff[i] * coeff[i];
        }
        int p_pos = 0;  // used pseudorange position
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats._sat[i_sat]._sys & SYS_ARRAY[i]) == SYS_ARRAY[i]) {
                p_pos = i; break;
            }
        if (_opt->_nsys == 1) p_pos = 0;
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

void CPntspp::GetAvdl(sat sats, double* sitepos, MatrixXd pos, MatrixXd &w) {
    int satnum = sats._nsats, num = 0;
    #pragma omp parallel for
    for (int i_sat = 0; i_sat < satnum; ++ i_sat) {
        if (!sats._sat[i_sat]._isused) continue;
        double dist = 0, coeff[3] = {0};
        for (int i = 0; i < 3; ++ i) {
            coeff[i] = sitepos[i] - sats._sat[i_sat]._pos[i];
            dist += coeff[i] * coeff[i];
        }
        dist = sqrt(dist);
        double rdot = 0;
        for (int i = 0; i < 3; ++ i) 
            rdot += -coeff[i] * sats._sat[i_sat]._vel[i] / dist;
        int p_pos = 0;
        for (int i = 0; i < MAXSYS; ++ i) 
            if ((sats._sat[i_sat]._sys & SYS_ARRAY[i]) == SYS_ARRAY[i])
                p_pos = i;
        if (_opt->_nsys == 1) p_pos = 0;
        double temp = sats._sat[i_sat]._clk[1];
        w(num++, 0) = -sats._sat[i_sat]._obs->_D[0] * VEL_LIGHT / GetFreq(sats._sat[i_sat]._sys, FREQTYPE_L1) - 
                      rdot + temp * VEL_LIGHT - pos(3 + p_pos, 0);
    }
}

void CPntspp::evaluateAVD() {
    MatrixXd var_vel = _optimizer->GetVelVar();
    MatrixXd Q = var_vel.block<0, 0>(3, 3);
    double sigma = _optimizer->GetVelInternalSigma();

    for (int i = 0; i < 3; ++ i)
        _res->_sigma_vel[i] = sqrt(Q(i, i));
}
