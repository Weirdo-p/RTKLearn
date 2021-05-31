#ifndef _LEASTSQ_H_
#define _LEASTSQ_H_
#include "navigation/navicommon.h"
#include "navigation/matrix.h"

class CLeastsq {
public:
    CLeastsq();

public:
    bool optimize(MatrixXd B, MatrixXd P, MatrixXd w);
    MatrixXd Getx();
    
private:
    MatrixXd x_;
    
};
#endif // _LEASTSQ_H_