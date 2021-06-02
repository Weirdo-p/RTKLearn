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
            memset(rnxopt_.obstypepos_, -1, sizeof(unsigned short) * MAXSYS * MAXFREQ * 4);
            int obsnum = scanobsfile(in);
            in.clear(); in.seekg(0);
            readobsh(in, opt);
            readobsb(in, obsnum, opt);
        } else if (line[20] == 'N') {
            int satnums = scannav(in, opt);
            eph_->msg_ = new nav_t [satnums];
            eph_->num = satnums;
            memset(eph_->msg_, 0, sizeof(nav));
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
        if (p = strchr(sys_, line[0]))
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
                if (!(p = strchr(sys_, line[0])))
                    continue;
                int sys = int(p - sys_);
                int sysflag = code2sys(line[0]);
                if (((opt.navsys_ & SYS_ARRAY[i]) == SYS_ARRAY[i]) && (sysflag == SYS_ARRAY[i]))
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
        if ((freqtypecode & opt.freqtype_) != freqtypecode)
            continue;

        // unknown observation type
        if (!(p = strchr(obstype_, freqtype[0])))
            continue;
        obstypepos = int(p - obstype_);
        // unknown track mode
        if (!(p = strchr(mode_, freqtype[2])))
            continue;
        if (rnxopt_.obstype_[sys][freqpos * 4 + obstypepos].size() == 0) {
            rnxopt_.obstype_[sys][freqpos * 4 + obstypepos] = freqtype;
            rnxopt_.obstypepos_[sys][freqpos * 4 + obstypepos] = position;
        }
        if (maxfreqtps < freqnums && position < freqnums && i == maxfreqtps - 1)
            getline(in, line);
    }
    return 1;
}

int CDecodeRnx304::readobsb(ifstream &in, int obsnum, prcopt opt) {
    static unsigned short sitenum = 0;
    int epoch = 0;
    if (obss_[sitenum].obs_) return -1;
    obss_[sitenum].obs_ = new obs_t[obsnum];
    obss_[sitenum].obsnum_ = obsnum;
    obss_[sitenum].rcv_ = sitenum;
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
               obss_[sitenum].obs_[epoch].time, 0);
    return satnum;
}

bool CDecodeRnx304::decodeobsr(ifstream &in, int sitenum, int satnum, prcopt opt, int &epoch) {
    string line;
    const char* p;
    Sattime time = obss_[sitenum].obs_[epoch].time;
    for (int i = 0; i < satnum; ++i) {
        getline(in, line);
        if(!(p = strchr(sys_, line[0])))
            continue;
        int syspos = int(p - sys_);
        if ((SYS_ARRAY[syspos] & opt.navsys_) != SYS_ARRAY[syspos])
            continue;
        for (int pos = 0; pos < MAXFREQ * 4; ++ pos) {
            int obspos = rnxopt_.obstypepos_[syspos][pos];
            if (obspos == -1)
                continue;
            if (!(p = strchr(freqcode_[syspos], rnxopt_.obstype_[syspos][pos][1])))
                continue;
            int freqpos = int(p - freqcode_[syspos]);
            if (!(p = strchr(obstype_, rnxopt_.obstype_[syspos][pos][0])))
                continue;
            int modepos = int(p - obstype_);
            string s_obs, s_lli;
            obss_[sitenum].obs_[epoch].sys = code2sys(line[0]);
            obss_[sitenum].obs_[epoch].sat = str2num<int>(line.substr(1, 2));
            try {
                s_obs = line.substr(3 + obspos * 16, 14);
                s_lli = line.substr(17 + obspos * 16, 1);
            } catch (...) {
                memset(&(obss_[sitenum].obs_[epoch]), 0, sizeof(obs));
                continue;
            }
            obss_[sitenum].obs_[epoch].lli[freqpos * 4 + modepos] = str2num<int>(s_lli);
            obss_[sitenum].obs_[epoch].time = time;
            switch (modepos) {
            case 0: // pseudorange
                obss_[sitenum].obs_[epoch].P[freqpos] = str2num<double>(s_obs);
                break;
            case 1: // carrier phase
                obss_[sitenum].obs_[epoch].L[freqpos] = str2num<double>(s_obs);
                break;
            case 2: // doppler
                obss_[sitenum].obs_[epoch].D[freqpos] = str2num<double>(s_obs);
                break;
            case 3: // signal strength
                obss_[sitenum].obs_[epoch].S[freqpos] = str2num<double>(s_obs);
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
        if ((sysflag & opt.navsys_) != sysflag) continue;
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
    eph_->msg_[ephnum].sys_ = SYS_GPS;
    string value = line.substr(1, 2);
    eph_->msg_[ephnum].prn_ = str2num<double>(value);
    value = line.substr(4, 4); year = str2num<int>(value);
    value = line.substr(9, 2); month = str2num<int>(value);
    value = line.substr(12, 2); day = str2num<int>(value);
    value = line.substr(15, 2); hour = str2num<int>(value);
    value = line.substr(18, 2); min = str2num<int>(value);
    value = line.substr(21, 2); sec = str2num<double>(value);
    Common2Gps(Commontime(year, month, day, hour, min, sec), eph_->msg_[ephnum].toc_, 0);
    value = line.substr(23, 19); eph_->msg_[ephnum].clkbias_ = str2num<double>(value);
    value = line.substr(42, 19); eph_->msg_[ephnum].clkdrift_ = str2num<double>(value);
    value = line.substr(61, 19); eph_->msg_[ephnum].clkdrate_ = str2num<double>(value);
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
    eph_->msg_[ephnum].Iode_ = d_ephvalue[0]; eph_->msg_[ephnum].Crs_ = d_ephvalue[1];
    eph_->msg_[ephnum].Deltan_ = d_ephvalue[2]; eph_->msg_[ephnum].M0_ = d_ephvalue[3];
    // orbit-2
    eph_->msg_[ephnum].Cuc_ = d_ephvalue[4]; eph_->msg_[ephnum].ecc_ = d_ephvalue[5];
    eph_->msg_[ephnum].Cus_ = d_ephvalue[6]; eph_->msg_[ephnum].sqrtA_ = d_ephvalue[7];
    // orbit-3
    eph_->msg_[ephnum].toe_ = Sattime(d_ephvalue[18], d_ephvalue[8]);
    eph_->msg_[ephnum].Cic_ = d_ephvalue[9]; eph_->msg_[ephnum].Omega0_ = d_ephvalue[10];
    eph_->msg_[ephnum].Cis_ = d_ephvalue[11];
    // orbit-4
    eph_->msg_[ephnum].I0_ = d_ephvalue[12]; eph_->msg_[ephnum].Crc_ = d_ephvalue[13];
    eph_->msg_[ephnum].Omega_ = d_ephvalue[14]; eph_->msg_[ephnum].Omega_dot_ = d_ephvalue[15];
    // orbit-5
    eph_->msg_[ephnum].Idot_ = d_ephvalue[16];
    // orbit-6
    eph_->msg_[ephnum].SV_ = d_ephvalue[20]; eph_->msg_[ephnum].SVHealth_ = d_ephvalue[21];
    eph_->msg_[ephnum].Tgd_[0] = d_ephvalue[22]; eph_->msg_[ephnum].Iodc_ = d_ephvalue[23];
    // orbit-7
    eph_->msg_[ephnum].Tof_ = d_ephvalue[25];
    return 1;
}

int CDecodeRnx304::readbdseph(ifstream &in, string &line, prcopt opt, int ephnum) {
    int year, month, day, hour, min;
    double sec;
    eph_->msg_[ephnum].sys_ = SYS_BDS;
    string value = line.substr(1, 2);
    eph_->msg_[ephnum].prn_ = str2num<double>(value);
    value = line.substr(4, 4); year = str2num<int>(value);
    value = line.substr(9, 2); month = str2num<int>(value);
    value = line.substr(12, 2); day = str2num<int>(value);
    value = line.substr(15, 2); hour = str2num<int>(value);
    value = line.substr(18, 2); min = str2num<int>(value);
    value = line.substr(21, 2); sec = str2num<double>(value);
    Common2Gps(Commontime(year, month, day, hour, min, sec), eph_->msg_[ephnum].toc_, 0);
    eph_->msg_[ephnum].toc_ = eph_->msg_[ephnum].toc_ + BDT2GPST;
    value = line.substr(23, 19); eph_->msg_[ephnum].clkbias_ = str2num<double>(value);
    value = line.substr(42, 19); eph_->msg_[ephnum].clkdrift_ = str2num<double>(value);
    value = line.substr(61, 19); eph_->msg_[ephnum].clkdrate_ = str2num<double>(value);
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
    eph_->msg_[ephnum].Iode_ = d_ephvalue[0]; eph_->msg_[ephnum].Crs_ = d_ephvalue[1];
    eph_->msg_[ephnum].Deltan_ = d_ephvalue[2]; eph_->msg_[ephnum].M0_ = d_ephvalue[3];
    // orbit-2
    eph_->msg_[ephnum].Cuc_ = d_ephvalue[4]; eph_->msg_[ephnum].ecc_ = d_ephvalue[5];
    eph_->msg_[ephnum].Cus_ = d_ephvalue[6]; eph_->msg_[ephnum].sqrtA_ = d_ephvalue[7];
    // orbit-3
    eph_->msg_[ephnum].toe_ = Sattime(d_ephvalue[18], d_ephvalue[8]);
    eph_->msg_[ephnum].Cic_ = d_ephvalue[9]; eph_->msg_[ephnum].Omega0_ = d_ephvalue[10];
    eph_->msg_[ephnum].Cis_ = d_ephvalue[11];
    BDST2GPST(eph_->msg_[ephnum].toe_, eph_->msg_[ephnum].toe_);
    // eph_->msg_[ephnum].toe_ = eph_->msg_[ephnum].toe_ + BDT2GPST;
    // orbit-4
    eph_->msg_[ephnum].I0_ = d_ephvalue[12]; eph_->msg_[ephnum].Crc_ = d_ephvalue[13];
    eph_->msg_[ephnum].Omega_ = d_ephvalue[14]; eph_->msg_[ephnum].Omega_dot_ = d_ephvalue[15];
    // orbit-5
    eph_->msg_[ephnum].Idot_ = d_ephvalue[16];
    // orbit-6
    eph_->msg_[ephnum].SV_ = d_ephvalue[20]; eph_->msg_[ephnum].SVHealth_ = d_ephvalue[21];
    eph_->msg_[ephnum].Tgd_[0] = d_ephvalue[22]; eph_->msg_[ephnum].Tgd_[1] = d_ephvalue[23];
    // orbit-7
    eph_->msg_[ephnum].Tof_ = d_ephvalue[25]; eph_->msg_[ephnum].Iodc_ = d_ephvalue[26];
    return 1;
}

int CDecodeRnx304::scannav(ifstream &in, prcopt opt) {
    string line;
    int count = 0;
    while(!in.eof()) {
        getline(in, line);
        int sysflag = code2sys(line[0]);
        if ((sysflag & opt.navsys_) == sysflag)
            count ++;
    }
    eph_->num;
    return count;
}
