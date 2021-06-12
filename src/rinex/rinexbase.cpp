#include "navigation/rinex/rinexbase.h"
#include <fstream>
#include <algorithm>

CRnxBase::CRnxBase() {
    this->obss_ = nullptr;
    this->_eph = new nav;
    this->_eph->_msg = nullptr;
    this->nsites_ = 0;
}

CRnxBase::~CRnxBase() {
    if (nsites_ != 0) {
        for (int i = 0; i < nsites_; ++i) 
            if (obss_[i]._obs) 
                delete[] obss_[i]._obs;
        // obss_->_obs = nullptr;
        delete[] obss_;
        // obss_ = nullptr;
    }

    if(_eph) {
        if (_eph->_msg) {
            delete[] _eph->_msg; //_eph._msg = nullptr;
            _eph->_num = 0;
        }
        delete _eph; _eph = nullptr;
    }
}

string CRnxBase::readver(char* infile) {
    ifstream in(infile);
    if (!in) {
        cout << "open file error" << endl;
        exit(-1);
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
    string ver;
    if (label == "RINEXVERSION/TYPE") 
        ver =  line.substr(5, 4);
    in.close();
    return ver;
}

void CRnxBase::setsites(int sitenum) {
    if (nsites_ == 0) {
        obss_ = new obs[sitenum];
        nsites_ = sitenum;
    }
}

obs* CRnxBase::GetObs() {
    return this->obss_;
}

nav* CRnxBase::GetEph() {
    return this->_eph;
}

obs::obs() {
    this->_obs == nullptr;
}
