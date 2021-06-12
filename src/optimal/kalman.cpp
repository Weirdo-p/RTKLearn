#include "navigation/optimal/kalman.h"

CKalman::CKalman() {
    _dim = 0;
}

void CKalman::setState(MatrixXd state) {
    _state = state; _dim = state.row();
}
void CKalman::setVarSys(MatrixXd varsys) {
    _var_sys = varsys;
}

void CKalman::setVarObs(MatrixXd var_obs) {
    _var_obs = var_obs;
}

void CKalman::setStateTrans(MatrixXd state_trans) {
    _state_trans = state_trans;
}

void CKalman::setVarState(MatrixXd var_state) {
    _var_state = var_state;
}

void CKalman::setObsMatrix(MatrixXd design) {
    _design = design;
}

MatrixXd CKalman::optimize(MatrixXd obs, MatrixXd h) {
    int issuccess = 0;
    MatrixXd state_time_predict = _state_trans * _state;
    MatrixXd var_state_time_predict = _state_trans * _var_state * _state_trans.transpose() + _var_sys;
    MatrixXd K = var_state_time_predict * _design.transpose() * (
        _design * var_state_time_predict * _design.transpose() + _var_obs
    ).inverse(issuccess);
    _state = state_time_predict + K * (obs - h);
    MatrixXd temp = K * _design;
    MatrixXd identity(temp.row(), temp.col()); identity.Identity();
    _var_state = (identity - temp) * var_state_time_predict;
    return _state;
}

CKalman::~CKalman() {
    _dim = 0;
}