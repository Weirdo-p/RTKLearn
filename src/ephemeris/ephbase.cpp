#include "navigation/ephemeris/ephbase.h"

CEphBase::CEphBase() {

}

bool CEphBase::satclk(Sattime time, sat_s &sat) {
    if (!sat._eph) return false;
    double t = time._Sow - sat._eph->_toc._Sow, ts = t;
    // relfix(time, sat);
    for (int i = 0; i < 5; ++i) 
        t = ts - (sat._eph->_clkbias + t * sat._eph->_clkdrift + t * t * sat._eph->_clkdrate);
    sat._clk[0] = sat._eph->_clkbias + t * sat._eph->_clkdrift + t * t * sat._eph->_clkdrate;
    sat._clk[1] = sat._eph->_clkdrift + 2 * sat._eph->_clkdrate * t;
    return true;
}

bool CEphBase::satvel(Sattime time, sat_s &sat) {
    double dt = 1E-3;
    sat_s sat_prev = sat;
    satpos(time, sat_prev);
    sat_s sat_now = sat;
    Sattime tmp = time + dt;
    satpos(tmp, sat_now);
    for (int i = 0; i < 3; ++i) 
        sat._vel[i] = (sat_now._pos[i] - sat_prev._pos[i]) / dt;
}

Ellipsoid CEphBase::GetElli() {
    return etype_;
}
