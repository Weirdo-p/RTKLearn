add_library(navigation SHARED
        config.cpp
        coors.cpp
        utils.cpp
        timemodule.cpp
        navicommon.cpp
        atmosphere.cpp
        rinex/rinex304.cpp
        rinex/rinexbase.cpp
        ephemeris/ephbase.cpp
        ephemeris/ephgps.cpp
        ephemeris/ephbds.cpp
        pnt/pntbase.cpp
        pnt/pntrtk.cpp
        pnt/pntspp.cpp
        optimal/leastsq.cpp)