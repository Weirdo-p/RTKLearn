#include "navigation/navicommon.h"
#include <string.h>
#include <stdio.h>
sat_s::sat_s() {
    _eph = nullptr; _obs = nullptr;
}

sat::sat() {
    _nsats = 0;
    memset(this->_sat, 0, sizeof(sat_s) * MAXOBS);
}

sat_s::~sat_s() {
    _eph = nullptr; _obs = nullptr;
}

prcopt::~prcopt() {
    memset(this, 0, sizeof(prcopt));
}

res_t::~res_t() {
    memset(this, 0, sizeof(res_t));
}
