#ifndef _PNTRTK_H_
#define _PNTRTK_H_

#include "navigation/pnt/pntspp.h"

class CPntrtk : public CPntbase {
public: // constructors
    CPntrtk();
    CPntrtk(prcopt opt);
    ~CPntrtk();

public:
    /**************************
     * process positioning
    **************************/
    int process() override;
    bool rtk(sat* sats_epoch);
    static int obsnumber(sat* sats);
    virtual void outsol(Sattime time, int nobs) override;
    virtual void evaluate() override;
    virtual int inputobs(sat* sats) override;

protected: // helpers
    int setobs(obs* obss, int &base_pos, int &rover_pos, sat* sats);
    int setnav(nav* navs, sat* sat);
    void getddobs(sat* sats, double* sitepos, int* refsats, int* sysobs, MatrixXd &w);
    int excludesats(sat &sat);
    void updateres(res_t res);

protected:
    CPntspp* _spprunner;
};



#endif // _PNTRTK_H_