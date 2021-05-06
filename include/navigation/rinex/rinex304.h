#ifndef _RINEX_H_
#define _RINEX_H_

#include "navigation/navicommon.h"
#include "navigation/rinex/rinexbase.h"

class CDecodeRnx304 : public CRnxBase{
public: // constructors
    CDecodeRnx304();
    virtual ~CDecodeRnx304();

public: // main functions
    /*********************************************
     * decode rinex files
     * @param   infiles [in]    rinex files
     * @param   opt     [in]    processing options
    **********************************************/
    int decode(char* infiles, prcopt opt);

protected: // observation reader
    int readobsh(ifstream &in, prcopt opt);
    int readobsb(ifstream &in, int obsnum, prcopt opt);
    int readfreqtype(ifstream &in, string line, int sysflag, prcopt opt);

protected:
    int readephh(ifstream &in, prcopt opt);
    int readephb(ifstream &in, prcopt opt);
    int readgpseph(ifstream &in, string &line, prcopt opt, int ephnum);
    int readbdseph(ifstream &in, string &line, prcopt opt, int ephnum);

    /***************************************
     * scan file to get total nav satellite
     * @param   in  [in]    file stream
     * @param   opt [in]    processing options
     * @return  total sats
    ***************************************/
    int scannav(ifstream &in, prcopt opt);

protected: // helper
    /***********************************************
     * transfer char to navigation system (SYS_???)
     * @param   code    sys code
     * @return  system code (SYS_???)
    ***********************************************/
    int code2sys(char code);

    /************************************************
     * transfer code to frequency type(FREQTYPE_???)
     * @param   code    freq type
     * @return  frequency type
    ************************************************/
    int code2freqnum(int code);

    /***********************************************
     * scan rinex ovservation files to acquire 
     * total observation epochs
     * @param   in  rinex observation file
     * @return  total observation epochs
    ***********************************************/
    int scanobsfile(ifstream &in);

    /*************************************************
     * transfer string to epoch
     * @param   sitenum [in]    which site 
     * @param   epoch   [in]    epoch nums
     * @param   line    [in]    a line start with ">"
     * @return  satellite nums
    *************************************************/
    int decodeEpoch(int sitenum, int &epoch, string line);

    /************************************************
     * decode observation record
     * @param   in      [in]    file stream
     * @param   sitenum [in]    site number
     * @param   satnum  [in]    satellite nums
     * @param   epoch   [in]    epoch num
     * @return  true if success
    ************************************************/
    bool decodeobsr(ifstream &in, int sitenum, int satnum, prcopt opt, int &epoch);
};

#endif