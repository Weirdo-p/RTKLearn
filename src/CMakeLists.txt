add_library(navigation SHARED
        config.cpp
        coors.cpp
        utils.cpp
        timemodule.cpp
        rinex/rinex304.cpp
        rinex/rinexbase.cpp)