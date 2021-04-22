#ifndef _RINEXBASE_H_
#define _RINEXBASE_H_

#include "navigation/navicommon.h"

class CRnxBase {
public:
    CRnxBase();
    virtual ~CRnxBase();

protected:
    int             nsites_;    //  site number
    rnxopt          rnxopt_;    //  rinex options
    obs*            obss_;      //  observations

protected:
    const char freqcode_[MAXSYS][MAXFREQ+1] = {
        FREQCODE_GPS,  // GPS
        FREQCODE_BDS   // BDS
    };

    // code, carrier phase, doppler, signal strength
    const char* obstype_ = OBSTYPE;

    const char* mode_ = TRACKMODE;

    const char* sys_ = SYSTEMS;
};

#endif  // _RINEXBASE_H_