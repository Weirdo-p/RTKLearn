/*---------------------------------------------
   coors.h
   create on 13 Apr 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _COORS_H_
#define _COORS_H_
#include "navigation/navicommon.h"
/***********************************
 * transform radians to degree
 * @param   rad     [In]    degree
 * @return deg
***********************************/
double Rad2Deg(double rad);

/***********************************
 * transform degree to radians
 * @param   deg     [In]    degree
 * @return rad
***********************************/
double Deg2Rad(double deg);

/***********************************************************
 * transfer ecef to blh (radians)
 * @param   xyz         [in]     coordinates in ecef
 * @param   ellipsoid   [in]     ellipsoid type
 * @param   blh         [in/out] coordinate in blh (radians)
 * @return  true if success
***********************************************************/
bool XYZ2BLH(const double* xyz, Ellipsoid ellipsoid, double* blh);

/********************************************************
 * transfer ecef to neu
 * @param   station  [in]  station coordinates in ecef
 * @param   obj      [in]  object coordinates in ecef
 * @param   type     [in]  ellipsoid type
 * @param   neu      [in/out] coordinate in neu
 * @return true if success
********************************************************/
bool XYZ2NEU(const double* station, const double* obj, Ellipsoid type, double *neu);

#endif // _COOORS_H_