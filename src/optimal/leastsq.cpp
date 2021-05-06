#include "navigation/optimal/leastsq.h"


bool CLeastsq::optimize(MatrixXd B, MatrixXd P, MatrixXd w) {
    int istrue = 0;
    MatrixXd inv = B.transpose() * P * B;
    inv = inv.inverse(istrue);
    x_ = inv * B.transpose() * P * w;
    return bool(istrue);
}

MatrixXd CLeastsq::Getx() {
    return x_;
}

CLeastsq::CLeastsq() {

}