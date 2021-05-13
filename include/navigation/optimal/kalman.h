#ifndef _KALMAN_H_
#define _KALMAN_H_

#include "navigation/matrix.h"
#include "navigation/navicommon.h"

class CKalman {
public:
    CKalman();
    
private:
    int dim_;
    MatrixXd state_;
    MatrixXd var_obs_;
    MatrixXd var_sys_;
    MatrixXd var_state_;
    MatrixXd state_trans_;
    MatrixXd design_;
};

#endif // _KALMAN_H_