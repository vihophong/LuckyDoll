#ifndef BuildAIDAEvents_H
#define BuildAIDAEvents_H

#include <vector>
#include "TTree.h"
#include "TClonesArray.h"

#include "AIDA.h"
#include "rawaida.h"
#include "AIDAUnpackerGz.h"
#include "rawaida.h"

class BuildAIDAEvents
{
public:
    BuildAIDAEvents();
    virtual ~BuildAIDAEvents();

    //! Initilize
    void Init(char* aidafile);

    //! Set verbose
    void SetVerbose(int verbose){fverbose = verbose;}

    //! Set Gz file input
    void SetGzStream(){fisgzstream = true;}

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

    //! Set fast discriminator scan window
    void SetDiscriminatorTimeWindow(unsigned long long disWin) {
        if (disWin>0){
            fwindowDisc=disWin;
            aidaunpkg->EnableFastDiscriminator();
        }
    }

    //! Set (AIDAUnpacker) max timestamp offset between exTS and inTS deviation hit by hit
    void SetAIDAMaxTSOffset(long long maxtsoffset) {aidaunpkg->SetMaxOffSet(maxtsoffset);}

    //! Set Fill data or not
    void SetFillData(bool fill_flag){fflag_filldata = fill_flag;}

    //! Set Fill data or not
    void SetPulserInStream(bool pulser_in_stream){fflag_pulser_in_stream = pulser_in_stream;}

    //! Set Sum Energy Cut
    void SetSumEXCut(Double_t sumexcut[]){for (Int_t i=0;i<NumDSSD;i++) fsumexcut[i] = sumexcut[i];}
    void SetSumEYCut(Double_t sumeycut[]){for (Int_t i=0;i<NumDSSD;i++) fsumeycut[i] = sumeycut[i];}

    //! Set Correlation energy cut
    void SetEnergyCorrCut(Double_t corrcut){fcorrcut = corrcut;}

    //! Get total number of events
    int GetADNblock(){return fADNblock;}
    //! Get Current events
    int GetCurrentADBlock(){return fADcurblock;}
    //! Get aida beta events
    int GetCurrentBetaEvent(){return fADBetaEntry;}
    //! Get aida pulse events
    int GetCurrentPulserEvent(){return fADPulserEntry;}
    //! Get aida ion events
    int GetCurrentIonEvent(){return fADIonEntry;}

    //! Is this event beta?
    int IsBETA(){return fisbeta;}

    AIDA* GetAIDABeta(){return flocalaidaBETA;}

    AIDA* GetAIDAIon(){return flocalaidaION;}

    //! Read calibration table to memory
    void ReadCalibTable();

    //!Get Next AIDA Event
    bool GetNextEvent();

    //! Get the ION tree
    TTree* GetIonTree(){return fmtrION;}
    //! Get the BETA tree
    TTree* GetBetaTree(){return fmtrBETA;}
    //! Get the PULSER tree
    TTree* GetPulserTree(){return fmtrPULSER;}

    void BookTree(TTree* treeIon, TTree *treeBeta, TTree *treePulser, Int_t bufsize=32000);

    //! Close event
    bool CloseIonEvent();
    bool CloseBetaEvent();

private:
    //!verbose
    int fverbose;
    //! store last event (sum of all)
    int flastevent;

    //! is gz file in stream
    bool fisgzstream;

    int fflag_pulser_in_stream;

    //! flag if the avaiable event is beta event
    bool fisbeta;
    //! set no filling data
    bool fflag_filldata;

    char* fmappingfile;
    char* fthresholdfile;
    char* fcalibfile;

    //!  output tree for ION
    TTree* fmtrION;
    //!  output tree for BETA
    TTree* fmtrBETA;
    //!  output tree for PULSER
    TTree* fmtrPULSER;

    unsigned int fADBetaEntry;
    unsigned int fADIonEntry;
    unsigned int fADPulserEntry;

    unsigned int fADcurblock;

    //!aida number of entries
    unsigned int fADentries;
    unsigned int fADNblock;

    //!  aida(if source is a tree) tree
    TTree* fADtr;

    //! input object for AIDA
    AIDAUnpacker* aidaunpkg;

    //! time stamps

    unsigned long long fADts;
    unsigned long long fADtsBETA;
    unsigned long long fADtsION;


    //! output objects

    AIDA* faida;

    AIDA* flocalaidaBETA;
    AIDA* flocalaidaION;

    //! stuff for the merger
    unsigned long long fcurrentionts;
    unsigned long long fcurrentbetats;
    bool fflag_addFirstIonHit;
    bool fflag_addFirstBetaHit;


    //! transient time
    bool fflag_trans;

    //! store fast time stamp
    unsigned long long flastfastts[NumFee][NumChFee];
    unsigned long long flastfasttsEXT[NumFee][NumChFee];

    unsigned long long flastADIonts;
    unsigned long long flastADBetats;

    //! AIDA calibration table
    Double_t dssd_cal[NumDSSD][NumStrXY][2];

    //! event window for ION event
    unsigned long long fwindowIon;
    //! event window for BETA event
    unsigned long long fwindowBeta;

    //! window for fast discriminator scan
    unsigned long long fwindowDisc;

    //! aida sleep time
    unsigned long long faidatranst;

    //! Add ION and BETA hits
    void AddAIDAIonHits(rawaida_info aidaraw);
    void AddAIDABetaHits(rawaida_info aidaraw);

    //! Check fast discriminator mask (new MIDAS update)
    bool check_channel_mask(Int_t dssdNo, Int_t stripNo, Int_t pattern,Int_t fee_num);

    //! Sum energy cut
    Double_t fsumexcut[NumDSSD];
    //! Sum energy cut
    Double_t fsumeycut[NumDSSD];

    //! Correlation cut
    Double_t fcorrcut;

    unsigned int ftemp;
    rawaida_info aidaraw;


};

#endif // BuildAIDAEvents_H
