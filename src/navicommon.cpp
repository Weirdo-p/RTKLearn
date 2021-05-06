#include "navigation/navicommon.h"
#include <string.h>
#include <stdio.h>
sat_s::sat_s() {
    eph_ = nullptr; obs_ = nullptr;
}

sat::sat() {
    nsats_ = 0;
    memset(this->sat_, 0, sizeof(sat_s) * MAXOBS);
}

sat_s::~sat_s() {
    eph_ = nullptr; obs_ = nullptr;
}

prcopt::~prcopt() {
    memset(this, 0, sizeof(prcopt));
}
