#include <iostream>
#include "navigation/navicommon.h"
#include "navigation/config.h"
#include "navigation/timemodule.h"
#include "navigation/rinex/rinex304.h"
#include "navigation/utils.h"
using namespace std;

int main(int argv, char** argc) {
    if (argv <= 1) {
        cout << "no file input, please check" << endl;
        return 1;
    }
    int* a = nullptr;
    if (a) cout << "yes" << endl;
    else   cout << "no" << endl;
    char path[256] = CONFPATH;
    CConfig config(path);
    CDecodeRnx304 rnxdecoder;
    char* files[1024];
    int n = 0;
    for (int i = 1; i < argv - 1; ++ i)
        files[n++] = argc[i];
    rnxdecoder.decode(files, n, config.GetConf());
}