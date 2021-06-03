#include <iostream>
#include "navigation/navicommon.h"
#include "navigation/config.h"
#include "navigation/timemodule.h"
#include "navigation/rinex/rinex304.h"
#include "navigation/utils.h"
#include "navigation/ephemeris/ephbase.h"
#include "navigation/ephemeris/ephgps.h"
#include "navigation/ephemeris/ephbds.h"
#include "navigation/pnt/pntbase.h"
#include "navigation/pnt/pntspp.h"
#include "navigation/pnt/pntrtk.h"
#include "navigation/optimal/rtkekf.h"
#include "navigation/optimal/optimal.h"
#include "navigation/ambiguity/lambda.h"


using namespace std;

int main(int argv, char** argc) {
    if (argv <= 1) {
        cout << "no file input, please check" << endl;
        return 1;
    }
    char path[1024] = CONFPATH;
    CConfig config(path);

    char* files[1024];
    int n = 0;
    for (int i = 1; i < argv; ++ i)
        files[n++] = argc[i];
    CPntbase* pnt;
    if(config.GetConf().mode_ == MODE_RTK) 
        pnt = new CPntrtk(config.GetConf());
    else if (config.GetConf().mode_ == MODE_SINGLE) 
        pnt = new CPntspp(config.GetConf());
    pnt->readRinex(files, n);
    pnt->process();
    delete pnt;

    return 0;
}