#include "navigation/optimal/kalman.h"

CKalman::CKalman() {
    dim_ = 0;
}

void CKalman::setState(MatrixXd state) {
    state_ = state; dim_ = state.row();
}
void CKalman::setVarSys(MatrixXd varsys) {
    var_sys_ = varsys;
}

void CKalman::setVarObs(MatrixXd var_obs) {
    var_obs_ = var_obs;
}

void CKalman::setStateTrans(MatrixXd state_trans) {
    state_trans_ = state_trans;
}

void CKalman::setVarState(MatrixXd var_state) {
    var_state_ = var_state;
}

void CKalman::setObsMatrix(MatrixXd design) {
    design_ = design;
}

MatrixXd CKalman::filter(MatrixXd obs, MatrixXd h) {
    int issuccess = 0;
    MatrixXd state_time_predict = state_trans_ * state_;
    MatrixXd var_state_time_predict = state_trans_ * var_state_ * state_trans_.transpose() + var_sys_;
    MatrixXd K = var_state_time_predict * design_.transpose() * (
        design_ * var_state_time_predict * design_.transpose() + var_obs_
    ).inverse(issuccess);
    state_ = state_time_predict + K * (obs - h);
    MatrixXd temp = K * design_;
    MatrixXd identity(temp.row(), temp.col()); identity.Identity();
    var_state_ = (identity - temp) * var_state_time_predict;
    return state_;
}

CKalman::~CKalman() {
    dim_ = 0;
}