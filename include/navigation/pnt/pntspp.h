#ifndef _PNTSPP_H_
#define _PNTSPP_H_
#include "navigation/pnt/pntbase.h"

class CPntspp : public CPntbase {
public:
    CPntspp();
    CPntspp(prcopt opt);
    ~CPntspp();

public:
    virtual int process() override;
    int spp(sat* sat);
    int spp_site(int isite, sat* sats);
    virtual void avd(int isite, sat* sats) override;
    virtual int excludesats(sat &sat) override;
    virtual int inputobs(sat* sats) override;
    virtual void outsol(Sattime time, int nobs) override;
    virtual void evaluate() override;
    virtual void evaluateAVD() override;
    
public:
    res_t* getRes();

protected:
    /*****************************************
     * fill in design matrix
     * @param   sats            satellite
     * @param   B   [in/out]    design matrix
     * @return  
    *****************************************/
    void GetDesign(sat sats, double* sitepos, MatrixXd &B);

    void Getl(int i_site, sat sats, double* sitepos, MatrixXd pos, MatrixXd &w);

    int setobs(obs* obss, int &rover_pos, sat* sats);

    int setnav(nav* navs, sat* sat);

    int GetAvdDesign(sat sats, double* sitepos, MatrixXd &B);

    void GetAvdl(sat sats, double* sitepos, MatrixXd pos, MatrixXd &w);
};


#endif // _PNTSPP_H_