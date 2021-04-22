#ifndef _RINEX_H_
#define _RINEX_H_

#include "navigation/navicommon.h"
#include "navigation/rinex/rinexbase.h"

class CDecodeRnx304 : public CRnxBase{
public: // constructors
    CDecodeRnx304();

public: // main functions
    /*********************************************
     * decode rinex files
     * @param   infiles [in]    rinex files
     * @param   n       [in]    number of files
     * @param   opt     [in]    processing options
    **********************************************/
    int decode(char** infiles, int n, Prcopt opt);

private: // observation reader
    int readobsh(ifstream &in, Prcopt opt);
    int readobsb(ifstream &in, int obsnum, Prcopt opt);
    int readfreqtype(ifstream &in, string line, int sysflag, Prcopt opt);

private: // helper
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
    int decodeEpoch(int sitenum, int epoch, string line);

    /************************************************
     * decode observation record
     * @param   in      [in]    file stream
     * @param   sitenum [in]    site number
     * @param   satnum  [in]    satellite nums
     * @param   epoch   [in]    epoch num
     * @return  true if success
    ************************************************/
    bool decodeobsr(ifstream &in, int sitenum, int satnum, int epoch);
};

#endif