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
    int excludesats(sat &sat);

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

protected:
    void bindRinex(char* path);
    void bindEph(int sys);

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

    /***************************************************
     * get frequencies
     * @param   sys         [in]    systems (SYS_???)
     * @param   freqflag    [in]    freqtype (FREQ_???)
     * @return  frequency
    ****************************************************/
    double GetFreq(int sys, int freqflag);

    /******************************
     * trop fix
     * @param   sats    satellites
     * @return
    ******************************/
    void tropfix(sat& sats, MatrixXd &w, double H);

    void setres(double* sitepos_ecef, double* sitepos_blh, Ellipsoid type, MatrixXd result);

    void selectPos(double** xzy, double** blh, int i_site);

    void getSysObs(sat sats, int* nobs_sys);

    double getSatUserPos(sat_s sats, double* sitepos);

    /************************************************
     * get weight by elevation 
     * 
     * sigma^2 = sigma0^2 * (1 + alpha * cos^2(E))
     * @param   elev    [in]    elevation(rad)
     * @param   simga0  [in]    prior sigma (m)
     * @param   alpha   [in]    coeff
     * @return  sigma
    ************************************************/
    double weightbyelev(double elev, double sigma0, double alpha);

public: // set function
    void setopt(prcopt opt);

protected:
    res_t*      res_;       /* result */
    prcopt*     opt_;       /* options */
    CRnxBase*   rnx_;       /* rinex reader */
    CEphBase*   satpos_;    /* satellite postion */
    CLeastsq    optimizer_; /* least square */
    CRtkekf     kf_;        /* kalman filter */
};

#endif // _PNTBASE_H_