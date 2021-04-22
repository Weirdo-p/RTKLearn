#ifndef _UTILS_H_
#define _UTILS_H_
#include <sstream>
#include "navigation/navicommon.h"

template <class T>
T str2num(string line) {
    stringstream buff;
    T num;
    buff << line; buff >> num;
    return num;
}

#endif // _UTILS_H_