#ifndef _EPHBASE_H_
#define _EPHBASE_H_
#include "navigation/navicommon.h"

class CEphBase {
public:
    CEphBase();

public:
    /*******************************************************
     * satellite position 
     * @param   time [in]      time
     * @param   sat [in/out]   satellite informations
     * @return  true if success
    ********************************************************/
    virtual bool satpos(Sattime time, sat_s &sat) { }

    /*************************************************
     * satellite velocity
     * @param   time [in]       time
     * @param   sat [in/out]    satellite information
     * @return true if success
    *************************************************/
    virtual bool satvel(Sattime time, sat_s &sat);

    /***************************************
     * satellite clock bias
     * @param   time    [in]        observation time
     * @param   sat     [in/out]    sat
     * @return  true if success
    ***************************************/
    virtual bool satclk(Sattime time, sat_s &sat);

    /***************************************
     * relative effective fix
     * @param   time    [in]        observation time
     * @param   sat     [in/out]    sat
     * @return  true if success
    ***************************************/
    virtual bool relfix(Sattime time, sat_s &sat) { }

    /*****************************************
     * earth rotation fix
     * @param   deltat          time of flight
     * @param   sat    [in/out] satellite information
     * @return
    *****************************************/
    virtual void earthRotateFix(double deltat, sat_s &sat) { }

public:
    Ellipsoid GetElli();

protected:
    Ellipsoid etype_;
};

#endif // _EPHBASE_H_