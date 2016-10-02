#ifndef BUILDEVENTS_H
#define BUILDEVENTS_H

#include <vector>
#include "TTree.h"
#include "TClonesArray.h"
#include "TArtEventInfo.hh"
#include "TArtBeam.hh"

#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "rawaida.h"
#include "AIDAUnpacker.h"
#include "rawaida.h"

/*!
  A container to keep track of the timestamps and corresponding detectors
*/
struct detector{
  //! timestamp of the detector hit
  unsigned long long int TS;
  //! ID for the detector, 0 BigRIPS, 1 AIDA, 2 BELEN, 3 Clover
  int ID;
};

class BuildEvents
{
public:
    BuildEvents();
    virtual ~BuildEvents();

    //! Initilize
    void Init(TTree* brtr, TTree* bltr, char* aidafile);

    //! Set verbose
    void SetVerbose(int verbose){fverbose = verbose;}

    //! Set file path containing a Mapping table
    void SetMappingFile(char* mapf){fmappingfile = mapf;}

    //! Set file path containing a Threshold table
    void SetThresholdFile(char* thrf){fthresholdfile = thrf;}

    //! Set file path containing a Calibration table
    void SetCalibFile(char* calf){fcalibfile = calf;ReadCalibTable();}

    //! Set event window ion
    void SetEventWindowION(unsigned long long windowIon) {fwindowIon=windowIon;}

    //! Set event window beta
    void SetEventWindowBETA(unsigned long long windowBeta) {fwindowBeta=windowBeta;}

    //! Set AIDA transient time window beta
    void SetAIDATransientTime(unsigned long long aidatranst) {faidatranst=aidatranst;}


    //! Get total number of events
    int GetADNblock(){return fADNblock;}
    //! Get Current events
    int GetCurrentADBlock(){return fADcurblock;}
    //! Get aida beta events
    int GetCurrentBetaEvent(){return fADBetaEntry;}
    //! Get aida beta events
    int GetCurrentIonEvent(){return fADIonEntry;}

    //! Read calibration table to memory
    void ReadCalibTable();

    //! Read bigrips entry
    bool ReadBigRIPS();

    //! Read aida entry
    bool ReadAIDA();

    //! Read belen entry
    bool ReadBELEN();

    //! Read clover entry
    bool ReadClover();

    //! Merge all
    bool Merge();

    //! Get the merged ION tree
    TTree* GetIonTree(){return fmtrION;}
    //! Get the merged BETA tree
    TTree* GetBetaTree(){return fmtrBETA;}

    //! Close event
    void CloseIonEvent();
    void CloseBetaEvent();

    bool GetNextEvent();

private:
    //!verbose
    int fverbose;
    //! store last event (sum of all)
    int flastevent;

    char* fmappingfile;
    char* fthresholdfile;
    char* fcalibfile;

    //! merged output tree for ION
    TTree* fmtrION;
    //! merged output tree for BETA
    TTree* fmtrBETA;

    //! bigrips, belen and aida current entry
    unsigned int fBRentry;
    unsigned int fBLentry;
    unsigned int fCLentry;

    unsigned int fADBetaEntry;
    unsigned int fADIonEntry;

    unsigned int fADcurblock;

    //! bigrips, belen and aida number of entries
    unsigned int fBRentries;
    unsigned int fBLentries;
    unsigned int fCLentries;
    unsigned int fADentries;
    unsigned int fADNblock;

    //! bigrips, bellen, clovers and aida(if source is a tree) entry
    TTree* fBRtr;
    TTree* fBLtr;
    TTree* fCLtr;
    TTree* fADtr;

    //! input object for Bigrips
    TClonesArray* fBRbeam;
    TClonesArray* fBReventinfo;

    //! input object reserve for BELEN and Clover (change later)
    TClonesArray* fBLCLeventinfo;
    TClonesArray* fBLCLhits;

    //! input object for AIDA
    AIDAUnpacker* aidaunpkg;

    //! time stamps
    unsigned long long fBRts;
    unsigned long long fBLts;
    unsigned long long fADts;
    unsigned long long fADtsBETA;
    unsigned long long fADtsION;


    //! output objects
    TClonesArray* fbeaminfo;
    TClonesArray* fbeam;
    AIDA* faida;
    BELEN* fbelen;
    Clover* fclover;

    AIDA* flocalaidaBETA;
    AIDA* flocalaidaION;

    //! stuff for the merger
    unsigned long long fcurrentionts;
    unsigned long long fcurrentbetats;

    //! transient time
    bool fflag_trans;
    //! last time stamp from each
    unsigned long long flastBRts;
    unsigned long long flastBLts;

    unsigned long long flastADIonts;
    unsigned long long flastADBetats;

    //! AIDA calibration table
    Double_t dssd_cal[NumDSSD][NumStrXY][2];

    //! event window for ION event
    unsigned long long fwindowIon;
    //! event window for BETA event
    unsigned long long fwindowBeta;

    //! aida sleep time
    unsigned long long faidatranst;

    //! list of detector hits
    vector<detector*> fdetectors;

    void AddAIDAIonHits(rawaida_info aidaraw);
    void AddAIDABetaHits(rawaida_info aidaraw);

};

#endif // BUILDEVENTS_H
