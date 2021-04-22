#include "navigation/rinex/rinexbase.h"

CRnxBase::CRnxBase() {
    this->obss_ = nullptr;
    this->nsites_ = 0;
}

CRnxBase::~CRnxBase() {
    if (nsites_ != 0) {
        for (int i = 0; i < nsites_; ++i) 
            if (obss_[i].obs_) {
                delete[] obss_[i].obs_;
                obss_[i].obs_ = nullptr;
            }
        delete[] obss_;
        obss_ = nullptr;
    }
}

obs::obs() {
    this->obs_ == nullptr;
}
