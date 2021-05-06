/*---------------------------------------------
   time.h time related module
   create on 12 Apr 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _TIMEMODULE_H_
#define _TIMEMODULE_H_

#include <iostream>
#include "navigation/navicommon.h"
using namespace std;

/**************************************************
 * function: to judge whether Common Time is legal
 * param:
 * @param   CT   Common Time
 * @return       flag 
 *      true     legal
 *      false   illegal      
**************************************************/
bool isLegal(const Commontime CT);

/**************************************************
 * function: to judge whether GPS Time is legal
 * param:
 * @param   GPST   GPS Time
 * @return         flag 
 *       true     legal
 *       false   illegal      
**************************************************/
bool isLegal(const Sattime GPST);

/**************************************************
 * function: to judge whether GPS Time is legal
 * param:
 * @param   MJD   GPS Time
 * @return        true if legal
**************************************************/
bool isLegal(const Mjdtime MJD);

/*****************************************
 * function: transform Common Time to MJD
 * param:
 * @param UTC  [in]     Universal Time
 * @param MJD  [out]    Modified Julian Day
 * @return              true if legal
*****************************************/
bool Common2Mjd(const Commontime UT, Mjdtime &MJD);

/***********************************
 * function:
 * transform MJD to Common Time
 * param:
 * @param MJD    [in]    Modified Julian Day
 * @param UT     [out]   Universal Time
 * @return               status code
***********************************/
bool Mjd2Common(const Mjdtime MJD, Commontime &UT);

/*******************************************
 * function:
 * transform MJD to GPST
 * param:
 * @param [in]    MJD  Modified Julian Day
 * @param [out]   GPST 
 * @return flag   status code
 *      true     success
 *      false  deadly error      
*******************************************/
bool Mjd2Gps(const Mjdtime MJD, Sattime &GPST);

/*******************************************
 * function:
 * transform MJD to GPST
 * param:
 * @param [in]    GPST  
 * @param [out]   MJD  Modified Julian Day
 * @return flag   status code
 *      true     success
 *      false  deadly error      
*******************************************/
bool Gps2Mjd(const Sattime GPST, Mjdtime &MJD);

/*******************************************
 * function:
 * transform MJD to GPST
 * param:
 * @param [in]    CT  Common Time
 * @param [out]   GPST 
 * @return flag   status code
 *      true     success
 *      false  deadly error      
*******************************************/
bool Common2Gps(const Commontime CT, Sattime &GPST, int leap);

/*******************************************
 * function:
 * transform MJD to Day Of Year
 * param:
 * @param  CT  [in] Common Time
 * @param  DOY [out]
 * @return flag   status code
*******************************************/
bool Common2Doy(const Commontime CT, unsigned short int &DOY);

/***************************************
 * function: to transform gpst to bdst
 * @param gpst [in]   GPS time
 * @param bdst [out]  BDS time
 * @return flag   status code
***************************************/
bool GPST2BDST(const Sattime gpst, Sattime &bdst);

/***************************************
 * function: to transform gpst to bdst
 * @param gpst [in]   GPS time
 * @param bdst [out]  BDS time
 * @return flag   status code
 *         true     success
 *         false  deadly error      
***************************************/
bool BDST2GPST(const Sattime bdst, Sattime &gpst);

/*************************************
 * transfer GPS time to UTC
 * @param   gpst  [in]  GPS time
 * @param   com   [out] UTC
 * @param   leap  [in]  leap seconds
*************************************/
bool GPST2Common(const Sattime gpst, Commontime& com, int leap);

/*********************************
 * satellite time minus
 * @param   t1 
 * @param   t2
 * @return seconds
*********************************/
double Sattimediff(const Sattime t1, const Sattime t2);

#endif // _TIME_H_