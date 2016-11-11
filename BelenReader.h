#ifndef BELENREADER_H
#define BELENREADER_H

#include <vector>
#include "TTree.h"
#include "TClonesArray.h"
#include "TreeData.h"
#include "BELEN.h"
#include "Clover.h"

#include <fstream>
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TROOT.h"
#include "TLatex.h"

class BelenReader
{
public:
    BelenReader();
    virtual ~BelenReader();

    //! set mapping files
    void SetMapping(char* mapf){fmappingfile = mapf;}
    //! se mapping (currently being called by Init)
    void GetMapping();
    //! initialization
    void Init(char* belenfile);
    //! book tree if we want to store files
    void BookTree(TTree* treeNeutron, TTree *treeGamma, TTree *treeAnc, Int_t bufsize=32000);
    //! close file reader
    void CloseReader(){finfile->Close();}

    //! get next event (any: gamma or neutron or ancillary)
    bool GetNextEvent();
    //! get next neutron event
    bool GetNextNeutronEvent();
     //! get next gamma event
    bool GetNextGammaEvent();
     //! get next ancillary event
    bool GetNextAncEvent();


    //! Get total number of events from the input file
    int GetNEvents(){return fnentries;}
    //! Get current processing entry
    int GetCurrentEvent(){return fcurentry;}
    //! Get current processing neutron entry
    int GetCurrentNeutronEvent(){return fBLNeuEntry;}
    //! Get current processing gamma entry
    int GetCurrentGammaEvent(){return fBLGamEntry;}
    //! Get current processing ancillary entry
    int GetCurrentAncEvent(){return fBLAncEntry;}

    //! neutron hits
    BELENHit* GetNeutron(){return flocalNeutron;}
    //! gamma hits
    CloverHit* GetGamma(){return flocalGamma;}
    //! ancillary hits
    BELENHit* GetAnc(){return flocalAnc;}

    //! Get event type, return 1 if neutron, 2 if gamma and 3 if ancillary!
    Double_t GetEnergy(){return fE;}
    ULong64_t GetTimestamp(){return fT;}
    UShort_t GetId(){return fId;}
    UShort_t GetType(){return ftype;}
    UShort_t GetIndex1(){return fIndex1;}
    UShort_t GetIndex2(){return fIndex2;}
    UShort_t GetInfoFlag(){return fInfoFlag;}
    std::string GetName(){return fName;}


    //! special function from MC simulation!
    void PertubateHe3(UShort_t He3Id);
    void PertubateClover(UShort_t CrystalId);

protected:
    //! mapping file
    char* fmappingfile;
    //! number of entries in belen file
    unsigned long long fnentries;
    //! current entry in belen file
    unsigned long long fcurentry;
    //! current neutron entry
    unsigned long long fBLNeuEntry;
    //! current gamma entry
    unsigned long long fBLGamEntry;
    //! current ancillary  entry
    unsigned long long fBLAncEntry;

    //! current time stamp of gamma
    unsigned long long fBLtsGamma;
    //! current time stamp of neutron
    unsigned long long fBLtsNeutron;
    //! current time stamp of ancillary
    unsigned long long fBLtsAnc;

    //! data structure
    Double_t        fE;
    ULong64_t       fT;
    UShort_t        fId;
    UShort_t        ftype;
    UShort_t        fIndex1;
    UShort_t        fIndex2;
    UShort_t        fInfoFlag;
    std::string          fName;
    Double_t        fposX;
    Double_t        fposY;
    Double_t        fposZ;

    Double_t fHe3Id2posX[MaxID];
    Double_t fHe3Id2posY[MaxID];
    Double_t fHe3Id2posZ[MaxID];
    Double_t fHe3Id2diameter[MaxID];

    Double_t fCrystalId2posX[MaxID];
    Double_t fCrystalId2posY[MaxID];
    Double_t fCrystalId2posZ[MaxID];

    //! neutron hits
    BELENHit* flocalNeutron;
    //! gamma hits
    CloverHit* flocalGamma;
    //! ancillary hits
    BELENHit* flocalAnc;

    //! input files
    TFile* finfile;

    //! a tree of input root files data
    BrikenTreeData* ftreedataNeuron;
    BrikenTreeData* ftreedataGamma;
    BrikenTreeData* ftreedataAnc;

    //! input tree
    TTree* ftree;
    //!  output tree for Neutron
    TTree* fmtrNeutron;
    //!  output tree for Gamma
    TTree* fmtrGamma;
    //!  output tree for Anc
    TTree* fmtrAnc;

    TRandom rr;

    bool fflag_filldata;

};

#endif // BELENREADER_H
