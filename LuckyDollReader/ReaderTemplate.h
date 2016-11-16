#include <iostream>
#include <iomanip>
#include <string>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "AIDA.h"
#include "TVectorD.h"
#include "TSystem.h"
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class AIDAHit+;
#pragma link C++ class AIDACluster+;
#pragma link C++ class AIDA+;
#endif

class Reader{
public:
    Reader(){
        gROOT->ProcessLine(".L AIDA.h++");
        aida = new AIDA;
    }
    virtual ~Reader(){}
    void InitBeta(TString infile){
        TFile* f = new TFile(infile);
        f->GetObject("beta",fChain);
        fChain->SetBranchAddress("aida",&aida);
    }
    AIDA* aida;
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    virtual Int_t    GetEntry(Long64_t entry){
        return fChain->GetEvent(entry);
    }
    virtual Long64_t    GetEntries(){
        return fChain->GetEntries();
    }
    //virtual void Loop(TString infile);
};
