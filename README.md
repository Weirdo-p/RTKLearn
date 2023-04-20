# Simple RTK
- Support dual frequency RTK of BDS and GPS
- Support both EKF and LS methods


usage: 
``` bash
mkdir build
cd build
cmake ..
make -j4
cd ../bin
./navi ${base obs file} ${rover obs file} ${nav file}
```

before run that command, modify process.conf and sites.conf in config folder to 
satisfy your requests.
