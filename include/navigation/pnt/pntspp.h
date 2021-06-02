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
    virtual int excludesats(sat &sat) override;

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

};


#endif // _PNTSPP_H_