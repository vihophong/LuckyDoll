#ifndef BUILDDECAY_H
#define BUILDDECAY_H

#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "Beam.h"
#include "DataStruct.h"
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH1F.h"
#include <vector>

typedef struct {
    double x,y;
    unsigned short z;
} aidaSimpleStruct;

class BuildDecay
{
public:
    BuildDecay();

    virtual ~BuildDecay();

    void SetInputFile(char* input){finput = input;}
    void Init();
    void ReadImplant();
    void ReadBeta();
    void ReadNeutron();
    void ReadF11Beam();
    void DoBuildDecay();

    void ReadImplant2();
    void ReadBeta2();
    void DoBuildDecay2();
    void DoBuildDecay3();
    void DoBuildDecay4();
    void BookTree(TTree* treeDecay);

    void SetMode(Int_t Mode){fMode=Mode;}

    void ShowMeYourFile(TFile* f){foutfile = f;}

    TH1F* GetHist1(){return fh1;}
protected:
    Int_t fMode;
    char* finput;

    //! file to be read
    TFile* finputFile;


    //! input tree
    TTree* ftrImp;
    TTree* ftrBeta;
    TTree* ftrNeutron;
    TTree* ftrF11Beam;

    //! number of enries in each tree
    Long64_t fnentriesImp;
    Long64_t fnentriesBeta;
    Long64_t fnentriesNeutron;
    Long64_t fnentriesF11Beam;

    //! data read from stream
    Implant* fimplant;
    Beta* fbeta;
    Neutron* fneutron;
    Ancillary* fF11Beam;



    //! Declare maps and their iterators (timestamp+entry)
    std::multimap < unsigned long long,std::pair< unsigned int, aidaSimpleStruct> > fimplantMap;
    std::multimap < unsigned long long, unsigned int > fbetaMap;
    std::multimap < unsigned long long, unsigned int> fneutronMap;
    std::multimap < unsigned long long, unsigned int> fbeamMap;

    std::multimap < unsigned long long,std::pair< unsigned int, aidaSimpleStruct> >::iterator fimplantMap_it;
    std::multimap < unsigned long long, unsigned int >::iterator fbetaMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fneutronMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fbeamMap_it;


    std::multimap < unsigned long long,std::pair< unsigned int, aidaSimpleStruct> > fbetaMap2;
    std::multimap < unsigned long long, unsigned int > fimplantMap2;
    std::multimap < unsigned long long,std::pair< unsigned int, aidaSimpleStruct> >::iterator fbetaMap_it2;
    std::multimap < unsigned long long, unsigned int >::iterator fimplantMap_it2;

    std::multimap < unsigned long long, unsigned int > fF11BeamMap;
    std::multimap < unsigned long long, unsigned int >::iterator fF11BeamMap_it;


     //! tree to be filled
     TFile* foutfile;
     TTree* ftreedecay;

     //! data struct to be filled
     double fdeltaxy;
     unsigned short fimplanti;     
     Implant* flocalimp;
     Beta* flocalbeta;
     Neutrons* flocalneutron;


     unsigned short fnF11Beam;
     Ancillary* flocalF11Beam;


     //! timewindows
     long long fBetaImplantTWup;
     long long fBetaImplantTWlow;
     long long fBetaNeutronTWup;
     long long fBetaNeutronTWlow;
     //! pixelation
     double fmaxdeltaxy;

     long long fBetaF11BeamTWup;
     long long fBetaF11BeamTWlow;
     long long fIonF11BeamTWup;
     long long fIonF11BeamTWlow;

     //! temp
     TH1F* fh1;

};

#endif // BUILDDECAY_H
