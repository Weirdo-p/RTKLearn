/*---------------------------------------------
   navicommon.h
   create on 12 Apr 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _NAVICOMMON_H_
#define _NAVICOMMON_H_
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include "navigation/matrix.h"

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
#define FREQ_B1I    1561.098e6  /* BDS B1 */
#define FREQ_B1C    1575.42e6   /* BDS B1C / B1A */
#define FREQ_B2A    1176.45e6   /* BDS B2A */    
#define FREQ_B2     1207.140e6  /* BDS B2 / B2b */
#define FREQ_B2AB   1191.795e6  /* BDS B2(B2a + B2b) */
#define FREQ_B3I    1268.52e6   /* BDS B3 / B3A */ 

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
#define FREQCODE_BDS    "267"       /* BDS frequency in rinex 3.04 */
#define TRACKMODE       "CWI"       /* supported track mode */
#define OBSTYPE         "CLDS"      /* observation type */
#define SYSTEMS         "GC"        /* supported systems */

// ephemeris type --------------------------------------------------
#define EPH_PREC    0       /* using precise ephemris */
#define EPH_BRDC    1       /* using broadcast ephemeris */

// time related ----------------------------------------------------
#define GPST2UTC   -18      /* difference between GPS and UTC */
#define BDT2GPST    14      /* difference between BDS and GPS */

// satellite define ------------------------------------------------
#define MAXGPSSATS  32      /* MAX gps satellites number */
#define MAXBDSSATS  64      /* MAX bds satellites number */
#define MAXSTAS     (MAXGPSSATS + MAXBDSSATS)
#define MAXOBS      64      /* MAX observations on an epoch */

// max supported sites ---------------------------------------------
#define MAXSITES    2       /* max used sites */

// solution and optput ---------------------------------------------
#define SOLTYPE_FLOAT   0   /* output float solutions */
#define SOLTYPE_FIX     1   /* output fix solutions */
#define PROC_LS         0   /* output Least Square solution */
#define PROC_KF         1   /* output Extended Kalman Filter solution (forward) */

// ambiguity related --------------------------------------------------
#define RATIO_THRES     3   /* ambiguity fix threshold */
#define FIX_SOLU        1
#define FLOAT_SOLU      2
#define SPP_SOLU        5        

// Array declaration -----------------------------------------------
const int SYS_ARRAY[] = {SYS_GPS, SYS_BDS};
const int FREQ_ARRAY[] = {FREQTYPE_L1, FREQTYPE_L2, FREQTYPE_L3};
const int GEO[] = {1, 2, 3, 4, 5, 59, 60, 61};

// typedef matrix
typedef Matrix<double, 3, 3> Matrix3d;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, Dynamic, Dynamic> MatrixXd;


// ellipsoid type --------------------------------------------------
enum EllipsoidType { CGCS2000, WGS84 };   /* support WGS84 CGCS2000 */

// structure defination --------------------------------------------
struct prcopt {             /* processing options */
   int _mode;               /* mode for processing (_mode???) */
   int _navsys;             /* systems to use */
   int _freqtype;           /* frequency to use */
   int _ephtype;            /* broadcast or precise eph to use */
   int _clktype;            /* clock type(broadcast or precise) */
   int _soltype;            /* solution type 0: float, 1: fix */
   int _proctype;           /* processing type 0: LS, 1: KF */
   unsigned short _freqnum; /* number of used frequency */
   unsigned short _nsys;    /* number of systems */
   unsigned short _sitenum; /* number of site */
   double _elecutoff;       /* elevation cutoff (in radians) */
   double _base[3] = {0};   /* priori coordinates of base */
   double _rover[3] = {0};  /* priori coordinates of rover */
   string _nbase;           /* name of base */
   string _nrover;          /* name of rover */

public:
    ~prcopt();
};

struct Ellipsoid {          /* ellipsoid type(default WGS84) */
public:
    double _a;              /* major semi axis */
    double _b;              /* minor semi axis */
    double _c;              /* helper param */
    double _e2;             /* square of first eccentricity */
    double _alpha ;         /* oblateness */
    double _eprime2;        /* second eccentricity */
    double _miu;            /* gravity const */
    double _rotation;       /* earth rotation */
    EllipsoidType _type;    /* ellipsoid type */

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
   int     _Week;
   double  _Sow;

public: // constructor
   Sattime();
   Sattime(int week, double sow);

public: // overload
   friend ostream & operator<<(ostream &out, const Sattime GPST);
   Sattime operator-(const Sattime &a) const;
   Sattime operator-(const double &a) const;
   Sattime operator+(const double &a) const;
   bool operator!=(const Sattime &a) const;
   double _2sec();
};

struct rnxopt {   /* rinex options */
    string _obstype[MAXSYS][MAXFREQ * 4]; /* observation type */
    short _obstypepos[MAXSYS][MAXFREQ * 4]; /* position of observation type */
};

struct obs_t {              /* observations for an epoch */
    Sattime _time;           /* observation time */
    int _sys, _sat;           /* system, satellite id */
    int _lli[MAXFREQ * 4];   /* loss of lock indicator */
    int _S[MAXFREQ];         /* signal strength */
    double _L[MAXFREQ];      /* carrier phase observations (cycle) */
    double _P[MAXFREQ];      /* pseudorange observations (m) */
    double _D[MAXFREQ];      /* doppler observations */
};

struct obs {                /* observations for total */
    int _obsnum;            /* total observation numbers */
    int _rcv;               /* reciever number */
    obs_t* _obs;            /* observations */

public: // constructors
    obs();
};

struct nav_t {  /* emphemeris for a satellite in an epoch */
    Sattime _sig;                               /* signal broadcast time */
    Sattime _toc;                               /* ephemeris time (toc) */
    Sattime _toe;                               /* time of ephemeris */
    double _clkbias, _clkdrift, _clkdrate;      /* clock parameters */
    double _Iode, _Crs, _Deltan, _M0;           /* orbit-1 in rinex 3.04 */
                                                /* IODE is AODE for bds */
    double _Cuc, _ecc, _Cus, _sqrtA;            /* orbit-2 in rinex 3.04 */
    double _Cic, _Omega0, _Cis;                 /* orbit-3 in rinex 3.04 */
    double _I0, _Crc, _Omega, _Omega_dot;       /* orbit-4 in rinex 3.04 */
    double _Idot;                               /* orbit-5 in rinex 3.04 */
    double _SV, _SVHealth, _Tgd[3], _Iodc;      /* orbit-6 in rinex 3.04 */
                                                /* tgd[0] -- gps/bds->tgd1(b1/b3) */
                                                /* tgd[1] -- bds->tgd1(b2/b3) */
                                                /* tgd[2] -- reserved */
    double _Tof;                                /* orbit-7 in rinex 3.04 */
    int _sys, _prn;                             /* systems, prn number */
};

struct nav {    /* navigation message */
    nav_t* _msg;     /* a record */
    int _num;        /* total number */
};

struct sat_s {    /* satellite information */
    int _sys;           /* navigation system */
    int _prn;           /* sat id */
    double _pos[3];     /* satellite position ecef */
    double _vel[3];     /* satellite velocity ecef */
    double _elev;       /* elevation angle */
    double _azi;        /* azimuth angle */
    double _clk[2];     /* clock bias / drift */
    double _gf;         /* geometry-free detection method */
    bool   _isused;     /* true if using this satellite */
    obs_t* _obs;        /* observations */
    nav_t* _eph;        /* ephemeris */
public:
    sat_s();
    ~sat_s();
};

struct sat {    /* satellite for an epoch */
    int _nsats;             /* number of satellites */
    sat_s _sat[MAXOBS];     /* satellite information */

public:
    sat();
};

struct res_t {  /* result for an epoch */
    int    _ratio;                  /* LAMBDA ratio value */
    unsigned int _ambi_flag;        /* 1: fix solution, 2: float solution, 5: spp */
    double _rdop;                   /* rdop */
    double _sigma_neu[3];           /* internal reliability */
    double _sigma_vel[3];           /* velocity internal reliability */
    double _vel[3];                 /* velocity */
    double _bpos_ecef[3];           /* base position */
    double _rpos_ecef[3];           /* rover position */
    double _bpos_blh[3];            /* base position */
    double _rpos_blh[3];            /* rover position */
    double _baseline[3];            /* baseline result (n/e/u) */
    double _enu[3];                 /* e/n/u respectively */
    double _recv_clk[MAXSITES];     /* reciever clock */
public:
    ~res_t();
};

template <class T>
T str2num(string line) {
    stringstream buff;
    T num;
    buff << line; buff >> num;
    return num;
}

#endif // _NAVICOMMON_H_