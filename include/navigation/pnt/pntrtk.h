#ifndef _PNTRTK_H_
#define _PNTRTK_H_

#include "navigation/pnt/pntspp.h"

class CPntrtk : public CPntbase {
public: // constructors
    CPntrtk();
    CPntrtk(prcopt opt);
    ~CPntrtk();

public:
    /**************************
     * process positioning
    **************************/
    int process() override;

protected:
    /************************************
     * select observations and ephemeris
     * @param   sats    both base and rover site
    ************************************/
    int inputobs(sat* sats);

protected: // helpers
    int setobs(obs* obss, int &base_pos, int &rover_pos, sat* sats);
    int setnav(nav* navs, sat* sat);
    nav_t* searchnav(Sattime time, int sys, int prn, nav* navs);

    /********************************************************
     * select reference satellite
     * @param   refsat  [in/out]    reference satellites' prn
     *                              0: GPS, 1: BDS
     * @param   sats    [in]        satellites
     * @return  true if success  
    ********************************************************/
    int chooseref(sat sats, double* refsat);

    /*******************************************************
     * choose reference satellite for a specific system
     * @param   sats    [in]        satellites
     * @param   sysflag [in]        system flag
     * @return  reference satellite
    *******************************************************/
    int chooseref(sat sats, int sysflag);

    /******************************************************
     * get dimension of design matrix
     * @param   nobs    [in]        number of observations
     * @param   row     [out]       rows
     * @param   col     [out]       columns
    ******************************************************/
    void getDesignDim(int nobs, int &row, int &col);

    /******************************************************
     * get design matrix
     * @param   sats    [in]    satellite
     * @param   nobs    [in]    number of observations
     * @param   sitepos [in]    rover postion
     * @param   refsat  [in]    reference satellites
     * @param   B       [out]   design matrix
     * @return
    ******************************************************/
    void getDesign(sat* sats, int nobs, double* sitepos, double* refsats, MatrixXd &B);

    /***************************************************
     * search reference satellite
     * @param   satts   [in]    satellites in an epoch
     * @param   sysflag [in]    system flag 
     * @param   prn     [in]    prn number
     * @return  satellite information
    ***************************************************/
    sat_s findRef(sat sats, int sysflag, int prn);

protected:
    CPntspp* spprunner_;
};


#endif // _PNTRTK_H_