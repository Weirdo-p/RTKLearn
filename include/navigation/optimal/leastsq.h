#ifndef _LEASTSQ_H_
#define _LEASTSQ_H_
#include "navigation/navicommon.h"
#include "navigation/matrix.h"
#include "navigation/optimal/optimal.h"

class CLeastsq : public COptimal {
public:
    CLeastsq();
    CLeastsq(prcopt* opt);
    
public:
    bool optimize(MatrixXd B, MatrixXd P, MatrixXd w, int nobs);
    void getl(sat* sats, double* sitepos, int* refsats, MatrixXd pos, MatrixXd &w);
    MatrixXd Getx();
    bool optimize(sat* sats_epoch, res_t &res);
    virtual MatrixXd GetVar() override;
    virtual double GetInternalSigma() override;
    virtual MatrixXd GetState() override;

private:
    MatrixXd x_;
    MatrixXd Q_;
    double sigma_; 
};

#endif // _LEASTSQ_H_