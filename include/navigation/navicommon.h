/*---------------------------------------------
   navicommon.h
   create on 12 Apr 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _NAVICOMMON_H_
#define _NAVICOMMON_H_
#include <math.h>
#include <iostream>

using namespace std;

// satellite related defination ------------------------------------
#define SYS_NONE    0x00    /* navigation system: NONE */
#define SYS_GPS     0x01    /* navigation system: GPS */
#define SYS_BDS     0x02    /* navigation system: BDS */
#define MAXSYS      2       /* max supported systems */

// refer to RINEX V3.05 --------------------------------------------
#define FREQ_L1     1575.42e6   /* GPS L1 */
#define FREQ_L2     1227.60e6   /* GPS L2 */
#define FREQ_L5     1176.45e6   /* GPS L5 */
#define FREQ_B1     1561.098e6  /* BDS B1 */
#define FREQ_B1C    1575.42e6   /* BDS B1C / B1A */
#define FREQ_B2A    1176.45e6   /* BDS B2A */    
#define FREQ_B2     1207.140e6  /* BDS B2 / B2b */
#define FREQ_B2AB   1191.795e6  /* BDS B2(B2a + B2b) */
#define FREQ_B3     1268.52e6   /* BDS B3 / B3A */ 

// satellites defination -------------------------------------------
#define MAXGPSSATS  32          /* maximum of GPS satellites */
#define MAXBDSSATS  64          /* maximum of BDS satellites */
#define SATS        (MAXGPSSATS + MAXBDSSATS)   /* total satellites */

// math defination -------------------------------------------------
#define PI          (atan(1) * 4)   /* pi 3.1415926.... */
#define VEL_LIGHT   299792458.0     /* speed of light */

// processing mode -------------------------------------------------
#define MODE_SINGLE  0
#define MODE_RTK     1

// frequency configuration
#define FREQTYPE_L1 0x01                /* frequency type: L1/B1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/B2 */
#define FREQTYPE_L3 0x04                /* frequency type: L5/L3 */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */
#define MAXFREQ     3
// default paths ---------------------------------------------------
#define CONFPATH  "./config/"         /* configuration path */
#define RESPATH   "./result/"         /* result save path */
#define NAME_CONF "process.conf"      /* name of processing configuraion file */
#define NAME_SITE "sites.conf"        /* name of sites information file */

// observation defination ------------------------------------------
#define FREQCODE_GPS    "125"       /* GPS frequency in rinex 3.04 */
#define FREQCODE_BDS    "276"       /* BDS frequency in rinex 3.04 */
#define TRACKMODE       "CWI"       /* supported track mode */
#define OBSTYPE         "CLDS"      /* observation type */
#define SYSTEMS         "GC"        /* supported systems */

// ephemeris type --------------------------------------------------
#define EPH_PREC    0
#define EPH_BRDC    1

// time related ----------------------------------------------------
#define GPST2UTC   -18
#define BDT2GPST    14

const int SYS_ARRAY[] = {SYS_GPS, SYS_BDS};
const int FREQ_ARRAY[] = {FREQTYPE_L1, FREQTYPE_L2, FREQTYPE_L3};

// ellipsoid type --------------------------------------------------
enum EllipsoidType {CGCS2000, WGS84};   /* support WGS84 CGCS2000 */

// structure defination --------------------------------------------
struct Prcopt {             /* processing options */
   int mode_;               /* mode for processing (MODE_???) */
   int navsys_;             /* systems to use */
   int freqtype_;           /* frequency to use */
   int ephtype_;            /* broadcast or precise eph to use */
   int clktype_;            /* clock type(broadcast or precise) */
   unsigned short freqnum_; /* number of used frequency */
   unsigned short nsys_;    /* number of systems */
   unsigned short sitenum_; /* number of site */
   double elecutoff_;       /* elevation cutoff (in radians) */
   double base_[3];         /* priori coordinates of base */
   double rover_[3];        /* priori coordinates of rover */
   string nbase_;           /* name of base */
   string nrover_;          /* name of rover */
};

struct Ellipsoid {          /* ellipsoid type(default WGS84) */
public:
    double a_;              /* major semi axis */
    double b_;              /* minor semi axis */
    double c_;              /* helper param */
    double e2_;             /* square of first eccentricity */
    double alpha_ ;         /* oblateness */
    double eprime2_;        /* second eccentricity */
    EllipsoidType type_;    /* ellipsoid type */

public:
    Ellipsoid();
    Ellipsoid(EllipsoidType type);
    Ellipsoid(double a, double b);

public:
    /**********************************
     * set ellipsoid params
     * @param   type    ellipsoid type
     * @return  true if success
    **********************************/
    bool SetEllipsoidParam(EllipsoidType type);

private:
    /**********************************
     * calculate ellipsoid params
     * @param   a   major semi axis
     * @param   b   minor semi axis
     * @return true if success
    **********************************/
    bool CalculateParam(double a, double b);
};

/*************************
 * common time structure
 * @param Year  ushort
 * @param Month ushort
 * @param Day   ushort
 * @param Hour  ushort
 * @param Min   ushort
 * @param Sec   double
*************************/
struct Commontime {
   unsigned short int  Year_, Month_, Day_, Hour_, Min_;
   double              Sec_;

public:
   Commontime(int year, int month, int day, int hour, int min, double sec);

   Commontime();

   friend ostream & operator<<(ostream &out, const Commontime UT);
};

/**************************
 * Julian Day structure
 * @param Day      int
 * @param FracDay  double
**************************/
struct Mjdtime {
   int     Day_;
   double  FracDay_;

   friend ostream & operator<<(ostream &out, const Mjdtime MJD);
};

/************************************
 * GPS Time structure
 * @param Week ushort
 * @param Sow  Second of Week double
************************************/
struct Sattime {
   int     Week_;
   double  Sow_;

   // default constructor
   Sattime();
   // overload operators
   friend ostream & operator<<(ostream &out, const Sattime GPST);
   Sattime operator-(const Sattime &a) const;
   Sattime operator-(const double &a) const;
   Sattime operator+(const double &a) const;
};

struct rnxopt {   /* rinex options */
    string obstype_[MAXSYS][MAXFREQ * 4]; /* observation type */
    short obstypepos_[MAXSYS][MAXFREQ * 4]; /* position of observation type */
};

struct obs_t {              /* observations for an epoch */
    Sattime time;           /* observation time */
    int sys, sat;           /* system, satellite id */
    int lli[MAXFREQ * 4];   /* loss of lock indicator */
    int S[MAXFREQ];         /* signal strength */
    double L[MAXFREQ];      /* carrier phase observations (cycle) */
    double P[MAXFREQ];      /* pseudorange observations (m) */
    double D[MAXFREQ];      /* doppler observations */
};

struct obs {                /* observations for total */
    int obsnum_;            /* total observation numbers */
    int rcv_;                /* reciever number */
    obs_t* obs_;            /* observations */

public:
    obs();
};

#endif // _NAVICOMMON_H_