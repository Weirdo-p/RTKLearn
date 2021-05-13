#include "navigation/coors.h"
#include "navigation/navicommon.h"
#include "navigation/matrix.h"

double Deg2Rad(double deg) {
    return (deg * PI) / 180.0;
}

double Rad2Deg(double rad) {
    return rad * 180.0 / PI;
}

Ellipsoid::Ellipsoid() {
    type_ = WGS84;
    SetEllipsoidParam(type_);
}

Ellipsoid::Ellipsoid(EllipsoidType type) {
    type_ = type;
    SetEllipsoidParam(type);
}

bool Ellipsoid::SetEllipsoidParam(EllipsoidType type) {
    switch (type) {
    case WGS84: {
        a_ = 6378137.0;
        b_ = 6356752.3142;
        miu_ = 3.986005e14;
        rotation_ = 7.2921151467e-5;
        CalculateParam(a_, b_);
        break;
    }
    case CGCS2000: {
        a_ = 6378137.0;
        b_ = 6356752.3141;
        miu_ = 3.986004418e14;
        rotation_ = 7.2921150e-5;
        CalculateParam(a_, b_);
        break;
    }
    default: break;
    }
}

bool Ellipsoid::CalculateParam(double a, double b) {
    if (a == 0 || b == 0) return false;
    double e, e_prime;
    c_ = a * a / b;
    alpha_ = (a - b) / a;
    e = sqrt(a * a - b * b) / a;
    e_prime = sqrt(a * a - b * b) / b;
    e2_ = e * e;
    eprime2_ = e_prime * e_prime;

    return true;
}

bool XYZ2NEU(const double* station, const double* obj, Ellipsoid type, double *neu) {
    double station_blh[3] = {0};
    XYZ2BLH(station, type, station_blh);
    double B = station_blh[0], L = station_blh[1];
    Matrix3d rotation;
    rotation(0, 0) = -sin(B) * cos(L); rotation(0, 1) = -sin(B) * sin(L); rotation(0, 2) = cos(B);
    rotation(1, 0) = -sin(L); rotation(1, 1) = cos(L); rotation(1, 2) = 0;
    rotation(2, 0) = cos(B) * cos(L); rotation(2, 1) = cos(B) * sin(L); rotation(2, 2) = sin(B);
    Vector3d vec;
    for (int i = 0; i < 3; ++i) 
        vec(i, 0) = obj[i] - station[i];
    Vector3d result = rotation * vec;
    for (int i = 0; i < 3; ++i)
        neu[i] = result(i, 0);
    return true;
}

bool XYZ2BLH(const double* xyz, Ellipsoid ellipsoid, double* blh) {
    blh[1] = atan2(xyz[1], xyz[0]);
    // 迭代初值  度为单位
    unsigned short int iteration = 0;
    // in deg
    double B0 = 1;
    while(iteration != 100) {
        double sinB = sin(Deg2Rad(B0));
        double W = sqrt(1 - ellipsoid.e2_ * sinB * sinB);
        double N = ellipsoid.a_ / W;
        double H = xyz[2] / sin(Deg2Rad(B0)) - N * (1 - ellipsoid.e2_);
        double up = xyz[2] + N * ellipsoid.e2_ * sinB;
        double down = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
        blh[0] = atan(up / down);
        blh[0] = Rad2Deg(blh[0]);
        blh[2] = H;
        double error = abs(blh[0] - B0);
        if(error < 1e-20)
            break;
        B0 = blh[0];
        iteration ++;
    }
    blh[0] = Deg2Rad(blh[0]);
    return true;
}
