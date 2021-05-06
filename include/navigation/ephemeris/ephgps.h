#ifndef _EPHGPS_H_
#define _EPHGPS_H_

#include "navigation/ephemeris/ephbase.h"
#include "navigation/navicommon.h"

class CEphGps : public CEphBase {
public:
    CEphGps();

public:
    /*******************************************************
     * satellite position 
     * @param   time [in]      time
     * @param   sat [in/out]   satellite informations
     * @return  true if success
    ********************************************************/
    virtual bool satpos(Sattime time, sat_s &sat) override;

    /***************************************
     * relative effective fix
     * @param   time    [in]        observation time
     * @param   sat     [in/out]    sat
     * @return  true if success
    ***************************************/
    virtual bool relfix(Sattime time, sat_s &sat) override;

    /*****************************************
     * earth rotation fix
     * @param   deltat  time of flight
     * @param   sat     satellite information
     * @return
    *****************************************/
    virtual void earthRotateFix(double deltat, sat_s& sat) override;

protected:
    // mean motion (rad / sec)
    double Calculaten0(nav_t nav);
    // eccentric anomaly
    double CalculateEk(double Mk, nav_t nav);
    // true anomaly
    double CalculateVk(double Ek, nav_t nav);
    // position in orbital plane
    void CalOrbPos(double rk, double uk, double* orbpos);
    // position in ecef
    void CalEcefPos(double* orbpos, double omegak, double ik, double* ecefpos);
    // ephemeris expiration check
    bool isExpired(double tk);
};

#endif // _EPHGPS_H_