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
    virtual bool optimize(MatrixXd B, MatrixXd P, MatrixXd w, int nobs) override;
    void getl(sat* sats, double* sitepos, int* refsats, MatrixXd pos, MatrixXd &w);
    MatrixXd Getx();
    bool optimize(sat* sats_epoch, res_t &res);
    virtual MatrixXd GetVar() override;
    virtual double GetInternalSigma() override;
    virtual MatrixXd GetState() override;
    virtual bool optimizeAVD(MatrixXd B, MatrixXd P, MatrixXd w, int nobs) override;
    virtual MatrixXd GetVel() override;
    virtual MatrixXd GetVelVar() override;
    virtual double GetVelInternalSigma() override;
    
private:
    MatrixXd _x;
    MatrixXd _vel;
    MatrixXd _Q;
    MatrixXd _vel_Q;
    double _sigma; 
    double _vel_sigma;
};

#endif // _LEASTSQ_H_