#include "ephemeris/ephbase.h"

CEphBase::CEphBase() {

}

bool CEphBase::satclk(Sattime time, sat_s &sat) {
    double t = time.Sow_ - sat.eph_->toc_.Sow_, ts = t;
    // relfix(time, sat);
    for (int i = 0; i < 5; ++i) 
        t = ts - (sat.eph_->clkbias_ + t * sat.eph_->clkdrift_ + t * t * sat.eph_->clkdrate_);
    sat.clk[0] = sat.eph_->clkbias_ + t * sat.eph_->clkdrift_ + t * t * sat.eph_->clkdrate_;
    sat.clk[1] = sat.eph_->clkdrift_ + 2 * sat.eph_->clkdrate_ * t;
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
        sat.vel_[i] = (sat_now.pos_[i] - sat_prev.pos_[i]) / dt;
}

Ellipsoid CEphBase::GetElli() {
    return etype_;
}
