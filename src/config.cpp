#include "navigation/config.h"
#include "navigation/coors.h"
#include "navigation/navicommon.h"
#include <string.h>
#include <string>
#include <sstream>
#include <algorithm>

CConfig::CConfig(char* path) {
    LoadConfig(path);
}

void CConfig::LoadConfig(char* path) {
    Loadprcopt(path);
    LoadSites(path);
}

void CConfig::reset() {
    _opt._navsys = SYS_NONE;
    _opt._mode = MODE_SINGLE;
    _opt._freqtype = FREQTYPE_L1;
    _opt._freqnum = 1;
    _opt._elecutoff = 0;
}

void CConfig::Loadprcopt(char* path) {
    char path_conf[256];
    strcpy(path_conf, path);
    strcat(path_conf, NAME_CONF);
    
    ifstream in(path_conf);
    if (!in) {
        cerr << "open configuration file error! " << endl;
        exit(-1);
    }
    while (!in.eof()) {
        string line;
        getline(in, line);
        ParseConfLine(line);
    }
}

void CConfig::ParseConfLine(string line) {
    int eqpos = 0, commentpos = 0;
    string label, value;
    for (int i = 0; i < line.size(); ++i) {
        if (i == 0 && line[i] == '#') return;
        if (line[i] == '=') eqpos = i;
        if (line[i] == '#') commentpos = i;
    }
    label = line.substr(0, eqpos);
    if (commentpos == 0) commentpos = line.size();
    value = line.substr(eqpos + 1, commentpos - eqpos - 1);
    auto itor = remove_if(label.begin(), label.end(), ::isspace);
    label.erase(itor, label.end());
    itor = remove_if(value.begin(), value.end(), ::isspace);
    value.erase(itor, value.end());
    if (!ParseLabel(label, value)) {
        cout << "FATAL: UNSUPPORTED OPTIONS" << endl;
        exit(-1);
    }
}

bool CConfig::ParseLabel(string label, string value) {
    bool issuccess = true;
    if (label == "navsys")          return SetSys(value);
    else if (label == "elecutoff")  return SetCutOff(value);
    else if (label == "mode")       return SetMode(value);
    else if (label == "eph")        return SetEphType(value);
    else if (label == "clk")        return SetClkType(value);
    else if (label == "freq")       return SetFreq(value);
    else if (label == "soltype")    return SetSoltype(value);
    else if (label == "filtertype") return SetProctype(value);
    else {
        cout << "unsupported options: " << label << endl;
        return false;
    }
}

bool CConfig::SetSoltype(string value) {
    transform(value.begin(), value.end(), value.begin(), ::tolower);
    if (value == "float")
        _opt._soltype = SOLTYPE_FLOAT;
    else if (value == "fix")
        _opt._soltype = SOLTYPE_FIX;
    else return false;

    return true;
}

bool CConfig::SetProctype(string value) {
    transform(value.begin(), value.end(), value.begin(), ::tolower);
    if (value == "kalman")
        _opt._proctype = PROC_KF;
    else if (value == "epoch") 
        _opt._proctype = PROC_LS;
    else return false;

    return true;
}

bool CConfig::SetSys(string value) {
    _opt._navsys = SYS_NONE; _opt._nsys = 0;
    for (int i = 0; i < value.size(); ++i) {
        if (value[i] == 'G') {
            _opt._navsys |= SYS_GPS; _opt._nsys ++;
        }
        if (value[i] == 'C') {
            _opt._navsys |= SYS_BDS; _opt._nsys ++;
        }
    }
    if (_opt._navsys == SYS_NONE) {
        cout << "no navigation system input" << endl;
        return false;
    }
    return true;
}

bool CConfig::SetCutOff(string value) {
    stringstream buff;
    double elev;
    buff << value; buff >> elev;
    if(elev < 1e-3) return false;
    _opt._elecutoff = Deg2Rad(elev);
    return true;
}

bool CConfig::SetMode(string value) {
    transform(value.begin(), value.end(), value.begin(), ::tolower);
    if (value == "spp")  _opt._mode = MODE_SINGLE;
    else if (value == "rtk") _opt._mode = MODE_RTK;
    else return false;

    return true;
}

bool CConfig::SetEphType(string value) {
    if (value == "brdc")    _opt._ephtype = EPH_BRDC;
    else if (value == "prec") _opt._ephtype = EPH_PREC;
    else return false;

    return true;
}

bool CConfig::SetClkType(string value) {
    if (value == "brdc")    _opt._clktype = EPH_BRDC;
    else if (value == "prec") _opt._clktype = EPH_PREC;
    else return false;

    return true;
}

bool CConfig::SetFreq(string value) {
    _opt._freqtype = 0; _opt._freqnum = 0;
    for(int i = 0; i < value.size(); ++i) {
        if (value[i] == '1'){
            _opt._freqtype |= FREQTYPE_L1;
            _opt._freqnum ++;
        } else if (value[i] == '2') {
            _opt._freqtype |= FREQTYPE_L2;
            _opt._freqnum ++;
        } else if (value[i] == '3') {
            _opt._freqtype |= FREQTYPE_L3;
            _opt._freqnum ++;
        } else if (value[i] == ',') continue;
        else if (value[i] == 'L') continue;
        else return false;
    }
    return true;
}

void CConfig::LoadSites(char* path) {
    strcat(path, NAME_SITE);

    ifstream in(path);
    if(!in) {
        cout << "can not find sites configurations " << endl;
        return ;
    }
    while (!in.eof()) {
        string line;
        getline(in, line);
        ParseSiteLine(line);
    }
    in.close();
}

void CConfig::ParseSiteLine(string line) {
    if (line[0] != ' ') return ;
    if (line[0] == '#') return ;
    static int sitenum = 0;
    stringstream buff;
    buff << line;
    if (sitenum == 0) {
        buff >> _opt._nbase;
        _opt._nbase = _opt._nbase.substr(0, 4);
        buff >> _opt._base[0] >> _opt._base[1] >> _opt._base[2];
    } else {
        buff >> _opt._nrover;
        _opt._nrover = _opt._nrover.substr(0, 4);
        buff >> _opt._rover[0] >> _opt._rover[1] >> _opt._rover[2];
    }
    sitenum++;
    _opt._sitenum = sitenum;
}

prcopt CConfig::GetConf() {
    return _opt;
}
