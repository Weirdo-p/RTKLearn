#include "navigation/coors.h"
#include "navigation/navicommon.h"

double Deg2Rad(double deg) {
    return (deg * PI) / 180.0;
}

double Rad2Deg(double rad) {
    return rad * 180.0 / PI;
}
