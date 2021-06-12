#ifndef _PNTBASE_H_
#define _PNTBASE_H_

#include "navigation/navicommon.h"
#include "navigation/ephemeris/ephbase.h"
#include "navigation/rinex/rinexbase.h"
#include "navigation/optimal/leastsq.h"
#include "navigation/optimal/rtkekf.h"

class CPntbase {
public: // constructors
    CPntbase();
    CPntbase(prcopt opt);
    virtual ~CPntbase();

public:
    /**************************************
     * read rinex files
     * @param   path    file paths
     * @param   n       file numbers
    **************************************/
    void readRinex(char** path, int n);

    /**************************
     * process positioning
    **************************/
    virtual int process() { }

    /**************************************
     * exclude satellite
     * @param   sat satellite
     * @return  used satellites number
    **************************************/
    virtual int excludesats(sat &sat) { }

    /*****************************************
     * earth rotation fix
     * @param   sat     satellite information
     * @return
    *****************************************/
    virtual void earthRotateFix(sat* sat);

    /*********************************************
     * calculate satellite azimuth and elevation
     * @param   site    site position
     * @param   sat     satellite
     * @return  
    *********************************************/
    void satazel(double* site, sat &sat);

    /**************************
     * evalutate position solutions
    **************************/
    virtual void evaluate() { }

    /********************************
     * evaluate velocity solutions
    ********************************/
    virtual void evaluateAVD() { }
    
    /*******************************************************
     * find observations in an epoch
     * @param   sats    [in/out]    observations in an epoch
     * @return false if EOF
    *******************************************************/
    virtual int inputobs(sat* sats) { }

    /*****************************************
     * absolute velocity determination
     * @param   isite   site number
     * @param   sat     satellite information
    *****************************************/
    virtual void avd(int isite, sat* sats) { }

protected:
    void bindRinex(char* path);
    void bindEph(int sys);
    nav_t* searchnav(Sattime time, int sys, int prn, nav* navs);
    /**************************************
     * calculate satellite position
     * @param   sats    satellites batch
     * @return  true if success
    **************************************/
    virtual bool satpos(sat* sats);

    /**************************************
     * calculate satellite clock bias
     * @param   sats    satellites batch
     * @return  true if success
    **************************************/
    virtual bool satclk(sat* sats);

    /********************************
     * satellite velocity
     * @param   satellite batch
     * @return  true if success
    ********************************/
    virtual bool satvel(sat* sats);

    /************************************
     * search pseudoranges
     * @param   P   pseudoranges
     * @return  pseudorange
    ************************************/
    virtual double searchpseu(double* P);

    /**************************************
     * relative effect fix
     * @param   sats    [in]    satellite
     * @return  true if success
    ***************************************/
    virtual bool relativeeffect(sat* sats);

    /********************************
     * check systems used
     * @param   sats    satellite
     * @return  system numbers
    ********************************/
    int usesys(sat &sats);

    /*****************************************
     * get iono-free combination observations
     * @param   sat     satellite
     * @param   w       matrix
     * @return  combination observations
    ******************************************/
    void GetLC(sat sat, MatrixXd &w);

    /*****************************************
     * get iono-free combination observations
     * @param   sat     satellite
     * @param   w       matrix
     * @return  combination observations
    ******************************************/
    void GetNonCombine(sat sat, MatrixXd &w);

    /******************************
     * trop fix
     * @param   sats    satellites
     * @return
    ******************************/
    void tropfix(sat& sats, MatrixXd &w, double H);

    void setres(double* sitepos_ecef, double* sitepos_blh, Ellipsoid type, MatrixXd result);

    void selectPos(double** xzy, double** blh, int i_site);

    void getSysObs(sat sats, int* nobs_sys);

    
public: // static function
    /************************************************
     * get weight by elevation 
     * 
     * sigma^2 = sigma0^2 * (1 + alpha * cos^2(E))
     * @param   elev    [in]    elevation(rad)
     * @param   simga0  [in]    prior sigma (m)
     * @param   alpha   [in]    coeff
     * @return  sigma
    ************************************************/
    static double weightbyelev(double elev, double sigma0, double alpha);

    /***************************************************
     * get frequencies
     * @param   sys         [in]    systems (SYS_???)
     * @param   freqflag    [in]    freqtype (FREQ_???)
     * @return  frequency
    ****************************************************/
    static double GetFreq(int sys, int freqflag);

    static double getSatUserPos(sat_s sats, double* sitepos);

    /************************
     * try to fix ambiguity
    ************************/
    void fixambi();

    virtual void outsol(Sattime time, int nobs) { }
public: // set function
    void setopt(prcopt opt);

protected:
    res_t*      _res;       /* result */
    prcopt*     _opt;       /* options */
    CRnxBase*   _rnx;       /* rinex reader */
    CEphBase*   _satpos;    /* satellite postion */
    COptimal*   _optimizer; /* least square */
};

#endif // _PNTBASE_H_