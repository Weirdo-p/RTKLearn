#include "navigation/timemodule.h"
#include <iostream>
#include <iomanip>
using namespace std;

ostream & operator<<(ostream &out, const Commontime UT) {

    out << UT.Year_ << "   " << UT.Month_ << "   "  << UT.Day_
        <<  "   " << UT.Hour_ << "   "  << UT.Min_ << "   " << UT.Sec_ << endl;
    return out;
}

ostream & operator<<(ostream &out, const Mjdtime MJD) {
    out << MJD.Day_ <<  "   " << MJD.FracDay_ << endl;
    return out;
}

ostream & operator<<(ostream &out, const Sattime GPST) {
    out << GPST.Week_ << "   " << setprecision(15) << GPST.Sow_ << endl;
    return out;
}

Sattime::Sattime(int week, double sow) {
    Week_ = week; Sow_ = sow;
}


Sattime::Sattime() {
    Week_ = 0; Sow_ = 0;
}

Sattime Sattime::operator-(const Sattime &a) const {      
    Sattime r;
    r.Week_ = this->Week_ - a.Week_; r.Sow_ = this->Sow_ - a.Sow_;
    if(r.Sow_ < 0) {
        r.Sow_ += 604800;
        r.Week_ -= 1;
    }
    return r;
}

double Sattime::_2sec() {
    return Week_ * 604800.0 + Sow_;
}

Sattime Sattime::operator-(const double &a) const {
    Sattime r;
    r.Week_ = this->Week_; r.Sow_ = this->Sow_ - a;
    if(r.Sow_ < 0) {
        r.Sow_ += 604800.0;
        r.Week_ -= 1;
    }
    return r;
}

Sattime Sattime::operator+(const double &a) const {
    Sattime r;
    r.Week_ = this->Week_;
    r.Sow_ = this->Sow_ + a;
    if(r.Sow_ >= 604800) {
        r.Sow_ -= 604800;
        r.Week_ += 1;
    }
    return r;
}

bool isLegal(const Commontime CT) {
    if(CT.Day_ > 31 || CT.Hour_ > 24 || CT.Min_ > 60 ||
       CT.Sec_ > 60 || CT.Day_ < 0 || CT.Hour_ < 0 || 
       CT.Min_ < 0 || CT.Sec_ < 0)
        return false;
    else
        return true;
}

bool isLegal(const Sattime GPST) {
    if(GPST.Sow_ < 0 || GPST.Week_ < 0)
        return false;
    else
        return true;
}

bool isLegal(const Mjdtime MJD) {
    if(MJD.Day_ < 0 || MJD.FracDay_ >= 1)
        return false;
    else
        return true;
}

bool Common2Mjd(const Commontime UT, Mjdtime &MJD, int leap) {
    bool flag;
    unsigned short int y, m;

    // 输入合法性检测
    flag = isLegal(UT);
    if(!flag)
        return flag;

    if(UT.Month_ <= 2) {
        y = UT.Year_ - 1;
        m = UT.Month_ + 12;
    } else {
        y = UT.Year_;
        m = UT.Month_;
    }
    // double JD = int(365.25 * y) + int(30.6001 * (m + 1)) + UT.Day_ + 1720981.5 + (UT.Hour_ + UT.Min_ / 60.0 + (UT.Sec_ + leap) / 3600.0) / 24.0;
    MJD.Day_ = int(365.25 * y) + int(30.6001 * (m + 1)) + UT.Day_ + 1720981;
    MJD.FracDay_ = (UT.Hour_ + UT.Min_ / 60.0 + (UT.Sec_ + leap) / 3600.0) / 24.0;
    MJD.Day_ -= 2400000;
    // FracDay 合法性检测
    if(MJD.FracDay_ >= 1) {
        MJD.Day_ += int(MJD.FracDay_);
        MJD.FracDay_ -= int(MJD.FracDay_);
    }
    flag = true;

    return flag;
}

bool Mjd2Common(const Mjdtime MJD_, Commontime &UT, int leap) {
    bool flag = isLegal(MJD_);
    if (!flag)
        return flag;
    Mjdtime MJD = MJD_;
    MJD.FracDay_ += leap / 86400.0;
    // helpers
    double JD = MJD.Day_ + MJD.FracDay_ + 2400000.5;
    int a = int(JD + 0.5);
    int b = 1537 + a;
    int c = int((b - 122.1) / 365.25);
    int d = int(365.25 * c);
    int e = int((b - d) / 30.6001);
    double FracD = JD + 0.5;  

    FracD -= int(FracD);
    UT.Day_ = b - d - int(30.6001 * e) + FracD;
    UT.Month_ = e - 1 - 12 * int(e / 14.0);
    UT.Year_ = c - 4715 - int((7 + UT.Month_) / 10.0);

    UT.Hour_ = int(FracD * 24.0);
    UT.Min_ = int((FracD * 24.0 - UT.Hour_) * 60);
    UT.Sec_ = ((FracD * 24.0 - UT.Hour_) * 60 - UT.Min_) * 60;
    flag = true;

    return flag;
}

bool Mjd2Gps(const Mjdtime MJD, Sattime &GPST) {
    bool flag = isLegal(MJD);
    if(!flag)
        return flag;

    GPST.Week_ = int((MJD.Day_ + MJD.FracDay_ - 44244) / 7.0);
    GPST.Sow_ = (MJD.Day_ - GPST.Week_ * 7.0 - 44244.0) * 86400.0;
    GPST.Sow_ +=   MJD.FracDay_ * 86400.0;
    // 合法性判断
    flag = isLegal(GPST);

    return flag;
}

bool Gps2Mjd(const Sattime GPST, Mjdtime &MJD_) {
    bool flag = isLegal(GPST);
    if (!flag)
        return flag;

    double MJD = 44244.0 + GPST.Week_ * 7.0 + GPST.Sow_ / 86400.0;
    MJD_.Day_ = int(MJD);
    MJD_.FracDay_ = MJD - MJD_.Day_;

    return true;
}

bool Common2Gps(const Commontime CT, Sattime &GPST, int leap) {
    Mjdtime MJD;
    bool flag = Common2Mjd(CT, MJD, leap);
    if(!flag)
        return flag;
    flag = Mjd2Gps(MJD, GPST);
    return flag;
}

bool Common2Doy(const Commontime CT, unsigned short int &DOY) {
    // 当前时间的儒略日
    Mjdtime MJD_Now;
    bool flag = Common2Mjd(CT, MJD_Now, 0);
    if(!flag)
        return flag;

    // 当前年第一天的儒略日
    Commontime NowYear;
    Mjdtime MJD_Beg;
    NowYear.Year_ = CT.Year_;
    NowYear.Month_ = 1;
    NowYear.Day_ = 1;
    NowYear.Hour_ = 0;
    NowYear.Min_ = 0;
    NowYear.Sec_ = 0;
    flag = Common2Mjd(NowYear, MJD_Beg, 0);
    if(!flag)
        return flag;
    
    DOY = int(MJD_Now.Day_ + MJD_Now.FracDay_ - (MJD_Beg.Day_ + MJD_Beg.FracDay_)) + 1;
    return flag;
}

bool GPST2BDST(const Sattime gpst, Sattime &bdst) {
    bdst.Week_ = gpst.Week_ - 1356;
    bdst.Sow_ = gpst.Sow_ - BDT2GPST;
    if(bdst.Sow_ < 0) {
        bdst.Sow_ += 604800;
        bdst.Week_ -= 1;
    }
    return true;
}

bool BDST2GPST(const Sattime bdst, Sattime &gpst) {
    gpst.Week_ = bdst.Week_ + 1356;
    gpst.Sow_ = bdst.Sow_ + 14;
    if(gpst.Sow_ >= 604800 ) {
        gpst.Sow_ -= 604800;
        gpst.Week_ += 1;
    }
    return true;
}

bool GPST2Common(const Sattime gpst, Commontime& com, int leap) {
    Mjdtime mjd;
    Gps2Mjd(gpst, mjd);
    Mjd2Common(mjd, com, leap);
    return true;
}

Commontime::Commontime(int year, int month, int day, int hour, int min, double sec) {
    Year_ = year; Month_ = month; Day_ = day; Hour_ = hour;
    Min_ = min; Sec_ = sec;
}

Commontime::Commontime() {
    Year_ = Month_ = Day_ = Hour_ = Min_ = Sec_ = 0;
}

double Sattimediff(const Sattime t1, const Sattime t2) {
    Sattime diff = t1 - t2;
    return diff.Week_ * 604800.0 + diff.Sow_;
}

bool Sattime::operator!=(const Sattime &a) const {
    if (abs(Sattimediff(a, *this)))
        return true;
    return false;
}
