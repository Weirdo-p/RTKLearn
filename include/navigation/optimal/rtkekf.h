#ifndef _RTKEKF_H_
#define _RTKEKF_H_

#include "navigation/optimal/kalman.h"
#include "navigation/navicommon.h"
#include "navigation/optimal/optimal.h"

#define VAR_INIT_STATE  (30 * 30)

// RTK status -------------------------------------------------------
#define RTK_INIT        0x00    /* EKF without initialize */
#define RTK_NORMAL      0x01    /* EKF processing normally */
#define RTK_OBSCHANGE   0x02    /* satellite changed */
#define RTK_REFCHANGE   0x04    /* reference satellites changed */
#define RTK_LEAP        0x08    /* cycle slip detected */
#define RTK_STATUS_NUM  5       /* number of status */

// number of observation change
#define MAX_OBSCHANGE   5       /* max observation change */

// noise defination ------------------------------------------------
#define INIT_STATE_VAR  15      /* init state variance */

const int STATUS_ARRAY[] { RTK_INIT, RTK_NORMAL, RTK_LEAP, RTK_OBSCHANGE, RTK_REFCHANGE };

class CRtkekf : public CKalman, public COptimal {
public:
    CRtkekf(prcopt* opt);
    CRtkekf();
    ~CRtkekf();

public:
    /*******************************************************************
     * execute extended kalman filter to solve RTK problems
     * @param   sats_epoch  satellites
     * @param   res         results
     * @param   obs         DD observations
     * @param   nobs        number of observations
     * @param   refsats     reference satellites
     * @param   sysobs      number of each observations for each epoch
     * @return  state vectors
    *******************************************************************/
    bool optimize(sat* sats_epoch, res_t &res);

    /***************************************************************
     * get DD prediction observations by using time-predicted state
     * @param   sats    [in]    satellites
     * @param   state   [in]    state vectors (time predicted)
     * @param   b_pos   [in]    base positions 
     * @param   refsats [in]    reference satellites
     * @param   w       [out]   DD predicted observations
    ***************************************************************/
    void getDDApprox(sat* sats, MatrixXd state, double* b_pos, int* refsats, MatrixXd &w);

public: // initialize
    /**********************************************************
     * initialize processing covariance matrix (dialog matrix)
     * @param   dim     dimension of state vector
     * @param   var     variance 
    ***********************************************************/
    void initvarsys(int dim, double var);

    /*******************************************************************
     * update EKF status by checking observations, reference satellites,
     * week leap detections
     * @param   obs         not used for now, reserved
     * @param   sats_epoch  satellites
     * @param   nobs        number of observations
     * @param   refsats     reference satellites
     * @param   res         SPP results (rover and base position)
    *******************************************************************/
    void updateStatus(MatrixXd obs, sat* sats_epoch, int nobs, int* refsats, res_t res, int* sysobs);

    /***********************************************************
     * initialize state vetors, position parameters are set to 
     * SPP results, and ambiguity parameters are initialized by 
     * DD pseudorange and DD phase observations
     * 
     * @param   sats    satellites
     * @param   b_pos   base ecef positions
     * @param   refsats reference satellites
     * @param   res     SPP results
    ***********************************************************/
    void initstate(sat* sats, double* b_pos, int* refsats, res_t res);

    /*********************************************
     * initialize variance of state vectors 
     * @param   dim     dimension of state vectors
     * @param   var     variance of states
    *********************************************/
    void initvarstate(int dim, double var);

    /*****************************************************************************
     * check whether observation changed (increase or decrease)
     * @param   obssats_current  [in]      satellites
     * @param   increase         [out]     satellites changed information
     *                                     each column represents for a specific system. 
     *                                     each row represents for a changed satellite' s
     *                                     prn corresponding to the system which is appeared
     *                                     at current epoch. 0 is default value. 
     * @param   decrease         [out]     each column represents for a specific system. 
     *                                     each row represents for a changed satellite' s
     *                                     prn corresponding to the system which is lost
     *                                     at current epoch. 0 is default value.
     * @return  true if number of satellite changed
    *****************************************************************************/
    bool checkObsChange(MatrixXd obssats_current, MatrixXd* increase, MatrixXd* decrease);

    /******************************************************************
     * check number of observation increases
     * @param   obssats_current [in]    current observed satellites
     * @param   increase        [out]   new satellite at current epoch
     * @return  true if increased
    ******************************************************************/
    bool checkObsIncrease(MatrixXd obssats_current, MatrixXd* increase);

    /****************************************************************
     * check whether number of observations decreases
     * @param   obssats_current [in]    current observed satellites
     * @param   decrease        [out]   satellites lost
     * @return  true if decreases
    ****************************************************************/
    bool checkObsDecrease(MatrixXd obssats_current, MatrixXd* decrease);

    /*************************************************************************
     * initialize ambiguties
     * @param   ref_sat_r   [in]        reference satellite of rover
     * @param   ref_sat_b   [in]        reference satellite of base
     * @param   sat_r       [in]        satellite of rover
     * @param   sat_b       [in]        satellite of base
     * @param   ar_pos      [in/out]    ambiguities position in state matrix
    *************************************************************************/
    void initambiguity(sat_s ref_sat_r, sat_s ref_sat_b, sat_s sat_r, sat_s sat_b, int nddobs, int &ar_pos);

    void fixNumIncrease(sat* sats, int* obssys, int nobs, MatrixXd increase);

    void fixNumDecrease(sat* sats, int* obssys, int nobs, MatrixXd increase);

    bool havesat(sat_s sat, MatrixXd changeObs);

    bool havesat(MatrixXd obssats, MatrixXd changeObs, int &pos);

    /*********************************************************************
     * cycle slip detection using geometry free (GF combination)
     * @param   sats_epoch  [in]    satellite
     * 
     * the state will be re-initialized if the corresponding satellite is 
     * slipped 
    *********************************************************************/
    void gfcycle(sat* sats_epoch);

    /********************************************
     * search gf combination of last epoch
     * @param   sat_last    satellites last epoch
     * @param   obj         target satellite
     * @return  gf combination last epoch
    ********************************************/
    double searchgf(sat sat_last, sat_s obj);

public: // set function
    /***************************************
     * set processing options
     * @param   opt     processing options
    ***************************************/
    void setopt(prcopt* opt);

public: // get function
    virtual MatrixXd GetState() override;
    virtual MatrixXd GetVar() override;


private:
    MatrixXd obssats(sat* sats_epoch, int nobs, int* refsats);
    void findrefpos(sat* sats_epoch, int* refsats);
    double dist(double* a, double* b);
    void statusFix(sat* sats_epoch, int nobs, int* refsats, int* obssys);
    void refchange(int* refsats, int* obssys);
    void initstate_s(sat* sats, int prn, int sysflag, int* ddobs);
    int obsnumber(sat* sats);
    void getDesignDim(sat sats, int nobs, int &row, int &col);
    void getddobs(sat* sats, double* sitepos, int* refsats, int* sysobs, res_t res, MatrixXd &w);
    void sortRefChange(int* new_ref_position, int* refsats, int* obssys);
    
private:
    int LeapSatsPos_[MAXOBS];   // position in design matrix
    int refsats_[MAXSYS];       // reference satellites at last epoch
    int sysobs_[MAXSYS];        // number of observations for each epoch
    int nobs_;                  // number of observations at current epoch
    int status_;                // RTK status at current epoch
    MatrixXd obssats_;          // observed satellites at last epoch
                                // column 0: system flag, column 1: prn
                                // column 2: ref sats flag
    res_t results_;             // kalman results   
};

#endif // _RTKEKF_H_

