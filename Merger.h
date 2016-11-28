#ifndef MERGER_H
#define MERGER_H

#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "Beam.h"
#include "DataStruct.h"
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH1F.h"


class Merger
{
public:
    Merger();

    virtual ~Merger();

    void SetAIDAFile(char* aidafile){finputAida = aidafile;}
    void SetBigripsFile(char* bigripsfile){finputBigrips = bigripsfile;}
    void SetBrikenFile(char* brikenfile){finputBriken = brikenfile;}
    void Init();
    void ReadBigrips();
    void ReadAIDA();
    void ReadBRIKEN();
    void DoMergeBeta();
    void DoMergeNeutron();
    void DoMergeAnc();
    void BookTree(TTree* treeImplant, TTree* treeBeta, TTree* treeNeutron, TTree* treeAnc);

    //void InitTemp();
    //void BookTreeTemp(TTree* treeAnc);
    void InitBRIKEN();
    void BookTreeBRIKEN(TTree* treeAnc);
    void ReadTemp();
    void DoMergeImp();

    TH1F* GetHist1(){return fh1;}
protected:
    char* finputAida;
    char* finputBigrips;
    char* finputBriken;

    //! file to be read
    TFile* fAidaFile;
    TFile* fBigripsFile;
    TFile* fBrikenFile;
    TFile* fF11AncFile;

    //! tree to be read
    TTree* ftrAIDABeta;
    TTree* ftrAIDAIon;
    TTree* ftrBigrips;
    TTree* ftrNeutron;
    TTree* ftrGamma;
    TTree* ftrAnc;
    TTree* ftrF11Anc;

    //! number of enries in each tree
    Long64_t fnentriesAIDABeta;
    Long64_t fnentriesAIDAIon;
    Long64_t fnentriesBigrips;
    Long64_t fnentriesNeutron;
    Long64_t fnentriesGamma;
    Long64_t fnentriesAnc;
    Long64_t fnentriesF11Anc;


    //! data read from stream

    AIDA* faidaBeta;
    AIDA* faidaIon;
    Beam* fbigrips;
    CloverHit* fclover;
    BELENHit* fneutron;
    BELENHit* fanc;
    Ancillary* ff11anc;


    //! Declare maps and their iterators (timestamp+entry)
    std::multimap < unsigned long long, unsigned int> fbigripsMap;
    std::multimap < unsigned long long, std::pair<unsigned int,unsigned short> > faidaBetaMap;
    std::multimap < unsigned long long, unsigned int> faidaIonMap;
    std::multimap < unsigned long long, unsigned int> fhe3Map;
    std::multimap < unsigned long long, unsigned int> fcloverMap;
    std::multimap < unsigned long long, unsigned int> fancMap;
    std::multimap < unsigned long long, unsigned int>::iterator fbigripsMap_it;
    std::multimap < unsigned long long, std::pair<unsigned int,unsigned short> >::iterator faidaBetaMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator faidaIonMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fhe3Map_it;
    std::multimap < unsigned long long, unsigned int>::iterator fcloverMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fancMap_it;


    //! for Nishimura san
    std::multimap < unsigned long long, unsigned int> ff11ancMap;
    std::multimap < unsigned long long, unsigned int> ff11ancMap_it;

    std::multimap < unsigned long long, unsigned int> fF11LMap;
    std::multimap < unsigned long long, unsigned int> fF11RMap;
    std::multimap < unsigned long long, unsigned int> fVetoTopMap;
    std::multimap < unsigned long long, unsigned int> fVetoBotMap;
    std::multimap < unsigned long long, unsigned int> fVetoDownMap;
    std::multimap < unsigned long long, unsigned int>::iterator fF11MapL_it;
    std::multimap < unsigned long long, unsigned int>::iterator fF11MapR_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoTopMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoBotMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoDownMap_it;

    std::multimap < unsigned long long, unsigned int> fdETopMap;
    std::multimap < unsigned long long, unsigned int> fdEBotMap;
    std::multimap < unsigned long long, unsigned int>::iterator fdETopMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fdEBotMap_it;


     //! tree to be filled
     TTree* ftreeImplant;
     TTree* ftreeBeta;
     TTree* ftreeNeutron;
     TTree* ftreeAncillary;

     //! data struct to be filled
     Implant* flocalimp;
     Beta* flocalbeta;
     Neutron* flocalneutron;
     Ancillary* flocalancillary;

     //! timewindows
     long long fIonPidTWup;
     long long fIonPidTWlow;
     long long fIonNeutronTWup;
     long long fIonNeutronTWlow;
     long long fIonGammaTWup;
     long long fIonGammaTWlow;
     long long fIonAncTWup;
     long long fIonAncTWlow;

     long long fNeuGammaTWup;
     long long fNeuGammaTWlow;
     long long fNeuAncTWup;
     long long fNeuAncTWlow;

     long long fBetaGammaTWup;
     long long fBetaGammaTWlow;
     long long fBetaAncTWup;
     long long fBetaAncTWlow;

     long long fF11LRTWup;
     long long fF11LRTWlow;
     long long fF11LVetoTWup;
     long long fF11LVetoTWlow;
     long long fF11LVetoDownup;
     long long fF11LVetoDownlow;

     long long fF11LdEtopWup;
     long long fF11LdEbotWup;
     long long fIondETWup;
     long long fIondETWlow;

     //! temp
     TH1F* fh1;

};

#endif // MERGER_H
