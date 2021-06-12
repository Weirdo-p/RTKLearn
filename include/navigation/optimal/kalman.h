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
    MatrixXd optimize(MatrixXd obs, MatrixXd h);

protected:
    int _dim;               // dimension of states
    MatrixXd _state;        // state vector
    MatrixXd _var_obs;      // variance of observations
    MatrixXd _var_sys;      // variance of systems
    MatrixXd _var_state;    // variance of states
    MatrixXd _state_trans;  // transition matrix
    MatrixXd _design;       // design matrix
};

#endif // _KALMAN_H_