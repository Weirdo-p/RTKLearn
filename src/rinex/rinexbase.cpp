#include "navigation/rinex/rinexbase.h"
#include <fstream>
#include <algorithm>

CRnxBase::CRnxBase() {
    this->obss_ = nullptr;
    this->eph_ = new nav;
    this->eph_->msg_ = nullptr;
    this->nsites_ = 0;
}

CRnxBase::~CRnxBase() {
    if (nsites_ != 0) {
        for (int i = 0; i < nsites_; ++i) 
            if (obss_[i].obs_) 
                delete[] obss_[i].obs_;
        // obss_->obs_ = nullptr;
        delete[] obss_;
        // obss_ = nullptr;
    }

    if(eph_) {
        if (eph_->msg_) {
            delete[] eph_->msg_; //eph_.msg_ = nullptr;
            eph_->num = 0;
        }
        delete eph_; eph_ = nullptr;
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
    return this->eph_;
}

obs::obs() {
    this->obs_ == nullptr;
}
