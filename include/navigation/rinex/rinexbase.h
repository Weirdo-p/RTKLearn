#ifndef _RINEXBASE_H_
#define _RINEXBASE_H_

#include "navigation/navicommon.h"

class CRnxBase {
public:
    CRnxBase();
    virtual ~CRnxBase();

public: // main function
    /*********************************************
     * decode rinex files
     * @param   infile [in]    rinex files
     * @param   opt     [in]    processing options
    **********************************************/
    virtual int decode(char* infiles, prcopt opt) { }

public: // set function
    /******************************************
     * init observations
     * @param   sitenum [in]    number of site
    ******************************************/
    void setsites(int sitenum);

public: // get functions
    virtual obs* GetObs();
    virtual nav* GetEph();

public: // tools
    /*********************************************
     * read rinex version
     * @param   infile  [in]    input file
     * @return  rinex version in string
    *********************************************/
    string readver(char* infile);

protected:
    int             nsites_;        //  site number
    rnxopt          rnxopt_;        //  rinex options
    obs*            obss_;          //  observations
    nav*            _eph;           //  ephemeris

protected:
    // code, carrier phase, doppler, signal strength
    const char* _obstype = OBSTYPE;
    const char* _mode = TRACKMODE;
    const char* _sys = SYSTEMS;
};

#endif  // _RINEXBASE_H_