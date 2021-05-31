#ifndef _KALMAN_H_
#define _KALMAN_H_

#include "navigation/matrix.h"
#include "navigation/navicommon.h"

class CKalman { // basic kalman filter structure (EKF)
public:
    CKalman();
    virtual ~CKalman();
    

public: // set function
    void setState(MatrixXd state);
    void setVarSys(MatrixXd var);
    void setVarObs(MatrixXd var);
    void setStateTrans(MatrixXd state_trans);
    void setVarState(MatrixXd var_state);
    void setObsMatrix(MatrixXd obs_matrix);

public:
    /**************************************************
     * extended kalman filter
     * @param   obs [in]    observations
     * @param   h   [in]    initial value at x(k, k-1)
     * @return  state at time k
    **************************************************/
    MatrixXd filter(MatrixXd obs, MatrixXd h);

protected:
    int dim_;               // dimension of states
    MatrixXd state_;        // state vector
    MatrixXd var_obs_;      // variance of observations
    MatrixXd var_sys_;      // variance of systems
    MatrixXd var_state_;    // variance of states
    MatrixXd state_trans_;  // transition matrix
    MatrixXd design_;       // design matrix
};

#endif // _KALMAN_H_