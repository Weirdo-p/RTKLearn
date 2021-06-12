#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <algorithm>
#include <chrono>
#include "navigation/rinex/rinex304.h"
#include "navigation/utils.h"
#include "navigation/timemodule.h"

using namespace std;

CDecodeRnx304::CDecodeRnx304() {

}
CDecodeRnx304::~CDecodeRnx304() {
    
}

int CDecodeRnx304::decode(char* infiles, prcopt opt) {
    ifstream in(infiles);
    if(!in) {
        cout << "open " << infiles << " error" << endl;
        return -1;
    }
    string line, label;
    getline(in, line);
    try {
        label = line.substr(60, 20);
    } catch (...) {
        label = line.substr(60, line.size() - 60);
    }
    auto itor = remove_if(label.begin(), label.end(), ::isspace);
    label.erase(itor, label.end());
    if (label == "RINEXVERSION/TYPE") {
        string s_ver = line.substr(5, 4);
        if (s_ver != "3.04") {
            cout << "rinex version is " << s_ver << " which is not supported" << endl;
            return -1;
        }
        if (line[20] == 'O') {
            memset(&rnxopt_, 0, sizeof(rnxopt));
            memset(rnxopt_._obstypepos, -1, sizeof(unsigned short) * MAXSYS * MAXFREQ * 4);
            int obsnum = scanobsfile(in);
            in.clear(); in.seekg(0);
            readobsh(in, opt);
            readobsb(in, obsnum, opt);
        } else if (line[20] == 'N') {
            int satnums = scannav(in, opt);
            _eph->_msg = new nav_t [satnums];
            _eph->_num = satnums;
            memset(_eph->_msg, 0, sizeof(nav));
            in.clear(); in.seekg(0);
            readephh(in, opt);
            readephb(in, opt);
        }
    }
    in.close();
}

int CDecodeRnx304::scanobsfile(ifstream &in) {
    int count = 0;
    string line;
    const char* p;
    while (!in.eof()) {
        getline(in ,line);
        if (p = strchr(_sys, line[0]))
            count += 1;
    }
    return count;
}

int CDecodeRnx304::readobsh(ifstream &in, prcopt opt) {
    const char* p;
    while(!in.eof()) {
        string line, label;
        getline(in, line);
        try {
            label = line.substr(60, 20);
        } catch (...) {
            label = line.substr(60, line.size() - 60);
        }
        auto itor = remove_if(label.begin(), label.end(), ::isspace);
        label.erase(itor, label.end());

        if (label == "SYS/#/OBSTYPES") {
            for (int i = 0; i < MAXSYS; ++i) {
                if (!(p = strchr(_sys, line[0])))
                    continue;
                int sys = int(p - _sys);
                int sysflag = code2sys(line[0]);
                if (((opt._navsys & SYS_ARRAY[i]) == SYS_ARRAY[i]) && (sysflag == SYS_ARRAY[i]))
                    readfreqtype(in, line, sys, opt);
            }
        }
        if (label == "ENDOFHEADER") 
            break;
    }
    return 1;
}

int CDecodeRnx304::readfreqtype(ifstream &in, string line, int sys, prcopt opt) {
    const int maxfreqtps = 13;
    int position = 0;
    string s_freqnum = line.substr(3, 3);
    stringstream buff(s_freqnum);
    double freqnums; buff >> freqnums;
    buff.clear();
    // if (sysflag == SYS_GPS) sys = 0;
    // else if (sysflag == SYS_BDS) sys = 1;
    for (int i = 0; i < maxfreqtps; i++, position++) {
        int f, freqpos = -1, obstypepos = -1;
        const char* p;
        string freqtype = line.substr(7 + 4 * i, 3);
        buff << freqtype[1]; buff >> f;
        // unknown freq
        if (!(p = strchr(freqcode_[sys], freqtype[1])))
            continue;
        freqpos = int(p - freqcode_[sys]);
        int freqtypecode = code2freqnum(freqpos);
        if ((freqtypecode & opt._freqtype) != freqtypecode)
            continue;

        // unknown observation type
        if (!(p = strchr(_obstype, freqtype[0])))
            continue;
        obstypepos = int(p - _obstype);
        // unknown track mode
        if (!(p = strchr(_mode, freqtype[2])))
            continue;
        if (rnxopt_._obstype[sys][freqpos * 4 + obstypepos].size() == 0) {
            rnxopt_._obstype[sys][freqpos * 4 + obstypepos] = freqtype;
            rnxopt_._obstypepos[sys][freqpos * 4 + obstypepos] = position;
        }
        if (maxfreqtps < freqnums && position < freqnums && i == maxfreqtps - 1)
            getline(in, line);
    }
    return 1;
}

int CDecodeRnx304::readobsb(ifstream &in, int obsnum, prcopt opt) {
    static unsigned short sitenum = 0;
    int epoch = 0;
    if (obss_[sitenum]._obs) return -1;
    obss_[sitenum]._obs = new obs_t[obsnum];
    obss_[sitenum]._obsnum = obsnum;
    obss_[sitenum]._rcv = sitenum;
    while (!in.eof()) {
        string line;
        getline(in, line);
        if (line[0] == '>') {
            int satnum = 0;
            if (!(satnum = decodeEpoch(sitenum, epoch, line)))
                continue;
            decodeobsr(in, sitenum, satnum, opt, epoch);
        }
    }
    sitenum += 1;
}

int CDecodeRnx304::decodeEpoch(int sitenum, int &epoch, string line) {
    
    int year, month, day, hour, min, satnum, flag;
    double sec;
    string s = line.substr(31, 1);
    flag = str2num<int>(s);
    if (flag != 0) return 0;
    try {
        s = line.substr(2, 4); year = str2num<int>(s);
        s = line.substr(7, 2); month = str2num<int>(s);
        s = line.substr(10, 2); day = str2num<int>(s);
        s = line.substr(13, 2); hour = str2num<int>(s);
        s = line.substr(16, 2); min = str2num<int>(s);
        s = line.substr(19, 11); sec = str2num<double>(s);
        s = line.substr(32, 3); satnum = str2num<int>(s);
    } catch (...) {
        return 0;
    }
    Common2Gps(Commontime(year, month, day, hour, min, sec),
               obss_[sitenum]._obs[epoch]._time, 0);
    return satnum;
}

bool CDecodeRnx304::decodeobsr(ifstream &in, int sitenum, int satnum, prcopt opt, int &epoch) {
    string line;
    const char* p;
    Sattime time = obss_[sitenum]._obs[epoch]._time;
    for (int i = 0; i < satnum; ++i) {
        getline(in, line);
        if(!(p = strchr(_sys, line[0])))
            continue;
        int syspos = int(p - _sys);
        if ((SYS_ARRAY[syspos] & opt._navsys) != SYS_ARRAY[syspos])
            continue;
        for (int pos = 0; pos < MAXFREQ * 4; ++ pos) {
            int obspos = rnxopt_._obstypepos[syspos][pos];
            if (obspos == -1)
                continue;
            if (!(p = strchr(freqcode_[syspos], rnxopt_._obstype[syspos][pos][1])))
                continue;
            int freqpos = int(p - freqcode_[syspos]);
            if (!(p = strchr(_obstype, rnxopt_._obstype[syspos][pos][0])))
                continue;
            int modepos = int(p - _obstype);
            string s_obs, s_lli;
            obss_[sitenum]._obs[epoch]._sys = code2sys(line[0]);
            obss_[sitenum]._obs[epoch]._sat = str2num<int>(line.substr(1, 2));
            try {
                s_obs = line.substr(3 + obspos * 16, 14);
                s_lli = line.substr(17 + obspos * 16, 1);
            } catch (...) {
                memset(&(obss_[sitenum]._obs[epoch]), 0, sizeof(obs));
                continue;
            }
            obss_[sitenum]._obs[epoch]._lli[freqpos * 4 + modepos] = str2num<int>(s_lli);
            obss_[sitenum]._obs[epoch]._time = time;
            switch (modepos) {
            case 0: // pseudorange
                obss_[sitenum]._obs[epoch]._P[freqpos] = str2num<double>(s_obs);
                break;
            case 1: // carrier phase
                obss_[sitenum]._obs[epoch]._L[freqpos] = str2num<double>(s_obs);
                break;
            case 2: // doppler
                obss_[sitenum]._obs[epoch]._D[freqpos] = str2num<double>(s_obs);
                break;
            case 3: // signal strength
                obss_[sitenum]._obs[epoch]._S[freqpos] = str2num<double>(s_obs);
                break;
            default:
                break;
            }
        }
        epoch++;
    }
    return true;
}

int CDecodeRnx304::code2sys(char code) {
    if (code == 'G') return SYS_GPS;
    else if (code == 'C') return SYS_BDS;
    else return -1;
}

int CDecodeRnx304::code2freqnum(int pos) {
    if (pos == 0)       return FREQTYPE_L1;
    else if (pos == 1)  return FREQTYPE_L2;
    else if (pos == 2)  return FREQTYPE_L3;
    else return -1;
}

int CDecodeRnx304::readephh(ifstream &in, prcopt opt) {
    string line, label;
    while (!in.eof()) {
        getline(in, line);
        try {
            label = line.substr(60, 20);
        } catch (...) {
            label = line.substr(60, line.size() - 60);
        }
        auto itor = remove_if(label.begin(), label.end(), ::isspace);
        label.erase(itor, label.end());
        if (label == "ENDOFHEADER") break;
    }
    return 1;
}


int CDecodeRnx304::readephb(ifstream &in, prcopt opt) {
    int ephnum = 0;
    string line;
    while(!in.eof()) {
        getline(in, line);
        int sysflag = code2sys(line[0]);
        if ((sysflag & opt._navsys) != sysflag) continue;
        if (sysflag == SYS_GPS)
            readgpseph(in, line, opt, ephnum);
        else if (sysflag == SYS_BDS) 
            readbdseph(in, line, opt, ephnum);
        ephnum ++;
    }
}

int CDecodeRnx304::readgpseph(ifstream &in, string& line, prcopt opt, int ephnum) {
    int year, month, day, hour, min;
    double sec;
    _eph->_msg[ephnum]._sys = SYS_GPS;
    string value = line.substr(1, 2);
    _eph->_msg[ephnum]._prn = str2num<double>(value);
    value = line.substr(4, 4); year = str2num<int>(value);
    value = line.substr(9, 2); month = str2num<int>(value);
    value = line.substr(12, 2); day = str2num<int>(value);
    value = line.substr(15, 2); hour = str2num<int>(value);
    value = line.substr(18, 2); min = str2num<int>(value);
    value = line.substr(21, 2); sec = str2num<double>(value);
    Common2Gps(Commontime(year, month, day, hour, min, sec), _eph->_msg[ephnum]._toc, 0);
    value = line.substr(23, 19); _eph->_msg[ephnum]._clkbias = str2num<double>(value);
    value = line.substr(42, 19); _eph->_msg[ephnum]._clkdrift = str2num<double>(value);
    value = line.substr(61, 19); _eph->_msg[ephnum]._clkdrate = str2num<double>(value);
    double d_ephvalue[28] = {0};
    for (int i = 0; i < 7; ++i) {
        getline(in, line);
        for (int j = 0; j < 4; ++j) {
            if (i == 6 && j >= 1) break;
            value = line.substr(4 + j * 19, 19);
            d_ephvalue[i * 4 + j] = str2num<double>(value);
        }
    }
    
    // orbit-1
    _eph->_msg[ephnum]._Iode = d_ephvalue[0]; _eph->_msg[ephnum]._Crs = d_ephvalue[1];
    _eph->_msg[ephnum]._Deltan = d_ephvalue[2]; _eph->_msg[ephnum]._M0 = d_ephvalue[3];
    // orbit-2
    _eph->_msg[ephnum]._Cuc = d_ephvalue[4]; _eph->_msg[ephnum]._ecc = d_ephvalue[5];
    _eph->_msg[ephnum]._Cus = d_ephvalue[6]; _eph->_msg[ephnum]._sqrtA = d_ephvalue[7];
    // orbit-3
    _eph->_msg[ephnum]._toe = Sattime(d_ephvalue[18], d_ephvalue[8]);
    _eph->_msg[ephnum]._Cic = d_ephvalue[9]; _eph->_msg[ephnum]._Omega0 = d_ephvalue[10];
    _eph->_msg[ephnum]._Cis = d_ephvalue[11];
    // orbit-4
    _eph->_msg[ephnum]._I0 = d_ephvalue[12]; _eph->_msg[ephnum]._Crc = d_ephvalue[13];
    _eph->_msg[ephnum]._Omega = d_ephvalue[14]; _eph->_msg[ephnum]._Omega_dot = d_ephvalue[15];
    // orbit-5
    _eph->_msg[ephnum]._Idot = d_ephvalue[16];
    // orbit-6
    _eph->_msg[ephnum]._SV = d_ephvalue[20]; _eph->_msg[ephnum]._SVHealth = d_ephvalue[21];
    _eph->_msg[ephnum]._Tgd[0] = d_ephvalue[22]; _eph->_msg[ephnum]._Iodc = d_ephvalue[23];
    // orbit-7
    _eph->_msg[ephnum]._Tof = d_ephvalue[25];
    return 1;
}

int CDecodeRnx304::readbdseph(ifstream &in, string &line, prcopt opt, int ephnum) {
    int year, month, day, hour, min;
    double sec;
    _eph->_msg[ephnum]._sys = SYS_BDS;
    string value = line.substr(1, 2);
    _eph->_msg[ephnum]._prn = str2num<double>(value);
    value = line.substr(4, 4); year = str2num<int>(value);
    value = line.substr(9, 2); month = str2num<int>(value);
    value = line.substr(12, 2); day = str2num<int>(value);
    value = line.substr(15, 2); hour = str2num<int>(value);
    value = line.substr(18, 2); min = str2num<int>(value);
    value = line.substr(21, 2); sec = str2num<double>(value);
    Common2Gps(Commontime(year, month, day, hour, min, sec), _eph->_msg[ephnum]._toc, 0);
    _eph->_msg[ephnum]._toc = _eph->_msg[ephnum]._toc + BDT2GPST;
    value = line.substr(23, 19); _eph->_msg[ephnum]._clkbias = str2num<double>(value);
    value = line.substr(42, 19); _eph->_msg[ephnum]._clkdrift = str2num<double>(value);
    value = line.substr(61, 19); _eph->_msg[ephnum]._clkdrate = str2num<double>(value);
    double d_ephvalue[28] = {0};
    for (int i = 0; i < 7; ++i) {
        getline(in, line);
        for (int j = 0; j < 4; ++j) {
            if (i == 6 && j >= 1) break;
            value = line.substr(4 + j * 19, 19);
            d_ephvalue[i * 4 + j] = str2num<double>(value);
        }
    }
    // orbit-1
    _eph->_msg[ephnum]._Iode = d_ephvalue[0]; _eph->_msg[ephnum]._Crs = d_ephvalue[1];
    _eph->_msg[ephnum]._Deltan = d_ephvalue[2]; _eph->_msg[ephnum]._M0 = d_ephvalue[3];
    // orbit-2
    _eph->_msg[ephnum]._Cuc = d_ephvalue[4]; _eph->_msg[ephnum]._ecc = d_ephvalue[5];
    _eph->_msg[ephnum]._Cus = d_ephvalue[6]; _eph->_msg[ephnum]._sqrtA = d_ephvalue[7];
    // orbit-3
    _eph->_msg[ephnum]._toe = Sattime(d_ephvalue[18], d_ephvalue[8]);
    _eph->_msg[ephnum]._Cic = d_ephvalue[9]; _eph->_msg[ephnum]._Omega0 = d_ephvalue[10];
    _eph->_msg[ephnum]._Cis = d_ephvalue[11];
    BDST2GPST(_eph->_msg[ephnum]._toe, _eph->_msg[ephnum]._toe);
    // _eph->_msg[ephnum]._toe = _eph->_msg[ephnum]._toe + BDT2GPST;
    // orbit-4
    _eph->_msg[ephnum]._I0 = d_ephvalue[12]; _eph->_msg[ephnum]._Crc = d_ephvalue[13];
    _eph->_msg[ephnum]._Omega = d_ephvalue[14]; _eph->_msg[ephnum]._Omega_dot = d_ephvalue[15];
    // orbit-5
    _eph->_msg[ephnum]._Idot = d_ephvalue[16];
    // orbit-6
    _eph->_msg[ephnum]._SV = d_ephvalue[20]; _eph->_msg[ephnum]._SVHealth = d_ephvalue[21];
    _eph->_msg[ephnum]._Tgd[0] = d_ephvalue[22]; _eph->_msg[ephnum]._Tgd[1] = d_ephvalue[23];
    // orbit-7
    _eph->_msg[ephnum]._Tof = d_ephvalue[25]; _eph->_msg[ephnum]._Iodc = d_ephvalue[26];
    return 1;
}

int CDecodeRnx304::scannav(ifstream &in, prcopt opt) {
    string line;
    int count = 0;
    while(!in.eof()) {
        getline(in, line);
        int sysflag = code2sys(line[0]);
        if ((sysflag & opt._navsys) == sysflag)
            count ++;
    }
    _eph->_num;
    return count;
}
