add_library(common SHARED
        config.cpp
        coors.cpp
        utils.cpp
        timemodule.cpp
        navicommon.cpp
        atmosphere.cpp)

add_library(pnt SHARED
        pnt/pntbase.cpp
        pnt/pntrtk.cpp
        pnt/pntspp.cpp
        ephemeris/ephbase.cpp
        ephemeris/ephgps.cpp
        ephemeris/ephbds.cpp
        rinex/rinex304.cpp
        rinex/rinexbase.cpp)

add_library(optimal SHARED
        optimal/leastsq.cpp
        optimal/kalman.cpp
        optimal/rtkekf.cpp)
