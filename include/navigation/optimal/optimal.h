#ifndef _OPTIMAL_H_
#define _OPTIMAL_H_

#include "navigation/navicommon.h"

class COptimal {
public:
    COptimal();
    COptimal(prcopt* opt);
public: // set function
    void setopt(prcopt opt);
    /***************************************
     * set processing options
     * @param   opt     processing options
    ***************************************/
    void setopt(prcopt* opt);
    virtual bool optimize(sat* sats_epoch, res_t &res) { }

protected:
    void getDesignDim(sat sats, int nobs, int &row, int &col);
    void getDesign(sat* sats, int nobs, double* sitepos, int* refsats, MatrixXd &B);
    void getSysObs(sat sats, int* nobs);
    int findSysPos(int sysflag, int* obs_sys);
    int findFreqPos(int sysflag, int* obs_sys, int &last_freq, double &freq);
    void getweight(sat* sats, int* refsats, int nobs, MatrixXd &P);
    void getsubweight(sat_s sat_, MatrixXd &subcov);
    MatrixXd getSingleDiffCov(MatrixXd cov);
    MatrixXd getDoubleDiffCov(MatrixXd cov_sd, int* refpos, int* obs_sys, int nobs);
    sat_s findRef(sat sats, int sysflag, int prn);
    int chooseref(sat sats, int* refsat);
    int chooseref(sat sats, int sysflag);
    
protected:
    prcopt* opt_;  
};
#endif // _OPTIMAL_H_