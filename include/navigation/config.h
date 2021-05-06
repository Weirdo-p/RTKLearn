/*---------------------------------------------
   config.h
   create on 12 Apr 2021 ZHUOXU WHU
---------------------------------------------*/

#ifndef _CONFIG_H_
#define _CONFIG_H_
#include <fstream>
#include <iostream>
#include "navigation/navicommon.h"
using namespace std;

class CConfig {
public: // constructors
    CConfig();
    CConfig(char* path);

public: // main functions
    /***************************************************
     * load user configurations
     * @param   [in]    path    configuration file path
     * @return
    ***************************************************/
    void LoadConfig(char* path);

    /**************************************************
     * set configurations to defaults
     * @return
    ***************************************************/
    void reset();

    /*************************************************
     * get configurations
     * @return configurations
    *************************************************/
    prcopt GetConf();

private: // processing configurations 
    void Loadprcopt(char* path);    /* load process.conf */
    void ParseConfLine(string line);    /* parse a line of process.conf */
    bool ParseLabel(string label, string value);  /* parse parameter */
    void SetSys(string value);
    bool SetCutOff(string value);   // elevation cut-off
    bool SetMode(string value);     // processing mode
    bool SetEphType(string value);
    bool SetClkType(string value);
    bool SetFreq(string value);

private: // sites configurations
    void LoadSites(char* path);
    void ParseSiteLine(string line);

private:
    prcopt opt_;
};

#endif // _CONFIG_H_