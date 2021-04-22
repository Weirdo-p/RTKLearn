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
    LoadPrcopt(path);
    LoadSites(path);
}

void CConfig::reset() {
    opt_.navsys_ = SYS_NONE;
    opt_.mode_ = MODE_SINGLE;
    opt_.freqtype_ = FREQTYPE_L1;
    opt_.freqnum_ = 1;
    opt_.elecutoff_ = 0;
}

void CConfig::LoadPrcopt(char* path) {
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
    ParseLabel(label, value);
}

bool CConfig::ParseLabel(string label, string value) {
    bool issuccess = true;
    if (label == "navsys") SetSys(value);
    else if (label == "elecutoff") issuccess = SetCutOff(value);
    else if (label == "mode")   issuccess = SetMode(value);
    else if (label == "eph")    issuccess = SetEphType(value);
    else if (label == "clk")    issuccess = SetClkType(value);
    else if (label == "freq")   issuccess = SetFreq(value);
}

void CConfig::SetSys(string value) {
    opt_.navsys_ = SYS_NONE; opt_.nsys_ = 0;
    for (int i = 0; i < value.size(); ++i) {
        if (value[i] == 'G') {
            opt_.navsys_ |= SYS_GPS; opt_.nsys_ ++;
        }
        if (value[i] == 'C') {
            opt_.navsys_ |= SYS_BDS; opt_.nsys_ ++;
        }
    }
}

bool CConfig::SetCutOff(string value) {
    stringstream buff;
    double elev;
    buff << value; buff >> elev;
    if(elev < 1e-3) return false;
    opt_.elecutoff_ = Deg2Rad(elev);
    return true;
}

bool CConfig::SetMode(string value) {
    transform(value.begin(), value.end(), value.begin(), ::tolower);
    if (value == "single")  opt_.mode_ = MODE_SINGLE;
    else if (value == "rtk") opt_.mode_ = MODE_RTK;
    else return false;

    return true;
}

bool CConfig::SetEphType(string value) {
    if (value == "brdc")    opt_.ephtype_ = EPH_BRDC;
    else if (value == "prec") opt_.ephtype_ = EPH_PREC;
    else return false;

    return true;
}

bool CConfig::SetClkType(string value) {
    if (value == "brdc")    opt_.clktype_ = EPH_BRDC;
    else if (value == "prec") opt_.clktype_ = EPH_PREC;
    else return false;

    return true;
}

bool CConfig::SetFreq(string value) {
    opt_.freqtype_ = 0; opt_.freqnum_ = 0;
    for(int i = 0; i < value.size(); ++i) {
        if (value[i] == '1'){
            opt_.freqtype_ |= FREQTYPE_L1;
            opt_.freqnum_ ++;
        } else if (value[i] == '2') {
            opt_.freqtype_ |= FREQTYPE_L2;
            opt_.freqnum_ ++;
        } else if (value[i] == '3') {
            opt_.freqtype_ |= FREQTYPE_L3;
            opt_.freqnum_ ++;
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
    static int sitenum = 0;
    stringstream buff;
    buff << line;
    if (sitenum == 0) {
        buff >> opt_.nbase_;
        opt_.nbase_ = opt_.nbase_.substr(0, 4);
        buff >> opt_.base_[0] >> opt_.base_[1] >> opt_.base_[2];
    } else {
        buff >> opt_.nrover_;
        opt_.nrover_ = opt_.nrover_.substr(0, 4);
        buff >> opt_.rover_[0] >> opt_.rover_[1] >> opt_.rover_[2];
    }
    sitenum++;
    opt_.sitenum_ = sitenum;
}

Prcopt CConfig::GetConf() {
    return opt_;
}
