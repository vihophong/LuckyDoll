#ifndef MERGER_H
#define MERGER_H

#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "Beam.h"
#include "DataStructNew.h"
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"

#include "fstream"


#define MaxNRI 1000


const Int_t kMaxGamma = 50;
const Int_t kMaxNeutron = 50;

typedef struct {
    ULong64_t evt;
    ULong64_t ts; 	 //timestamp in ns
    Double_t t,x,y,ex,ey,ion_x,ion_y,ion_ex,ion_ey,zet,aoq,deltaxy;//ts diff in ns
    Short_t z,ion_z,multx,multy,multz,ndecay,nbeta;

    Int_t gc_hit;
    Double_t gc_E[kMaxGamma];
    Double_t gc_T[kMaxGamma];//gamma time in ns
    Int_t gc_ch[kMaxGamma];

    Int_t neu_hit;
    Double_t neu_E[kMaxNeutron];
    Double_t neu_T[kMaxNeutron];//moderation time in us
    Int_t neu_ch[kMaxNeutron];
    Double_t neu_x[kMaxNeutron];
    Double_t neu_y[kMaxNeutron];
    Int_t neub_hit;


} datatype;

typedef struct{
    int correntrybrips;
    int correntryf11r;
    int correntrydEtop;
    int correntrydEbot;
    int correntryvetodown;
}ImplantCorrelationVector;

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
    void ReadAIDA(unsigned int start = 0, unsigned int stop = 0);
    void ReadBRIKEN(unsigned int startN=0, unsigned int stopN=0,unsigned int startG=0, unsigned int stopG=0,unsigned int startA=0, unsigned int stopA=0);
    void DoMergeTClone();
    void DoMergeSingle();
    void DoMergeTest();
    void DoMergeYOnly();
    void BookTreeTClone(TTree* tree, TTree* treemlh, TTree* treemlhp1n, TTree* treemlhp2n, TTree *treemlhp1nb, TTree *treemlhp2nb);
    void BookTreeSingle(TTree* tree);
    void BookTreeNeutron(TTree* tree);
    void BookTreeImplant(TTree* tree);

    void SetNeutronOffsetTime(long long offset){fNeuBetaoffset=offset;}
    long long GetNeutronOffsetTime(){return fNeuBetaoffset;}
    void SetSumERankCut(Int_t rankcut){sumexyrankcut = rankcut;}
    Int_t GetSumERankCut(){return sumexyrankcut;}


    //! for pid separation
    void ReadPID(char* pidfile,Int_t ncutpts=20);
    void SetMergedFile(char* mergedfile){finputMerged = mergedfile;}
    void InitPIDSep();
    void BookPIDSepTree();
    void BookPIDSepSimpleTree();
    Int_t GetNri(){return nri;}
    TTree* GetTreeRI(Int_t i){if (i<0) return ftree; else return ftreeRI[i];}
    TTree* GetTreeImpRI(Int_t i){if (i<0) return ftreeimplantAll;return ftreeimplantRI[i];}

    TCutG* GetCUTRI(Int_t i){return cutg[i];}
    void DoSeparatePID();

    void DoSeparatePIDFinalTree();

    TH1F* GetHist1(){return fh1;}

    Long64_t GetF11VetoDeadtime(){return ff11vetodeadtime;}
    Long64_t GetDownstreamVetoDeadtime(){return fdownstreamvetodeadtime;}
    Long64_t GetFinalVetoDeadtime(){return fvetodeadtime;}

    Long64_t GetF11VetoTotaltime(){return ff11vetototaltime;}
    Long64_t GetDownstreamVetoTotaltime(){return fdownstreamvetototaltime;}
    Long64_t GetFinalVetoTotaltime(){return fvetototaltime;}

    void ResetSimpleData();

protected:
    char* finputAida;
    char* finputBigrips;
    char* finputBriken;

    //! file to be read
    TFile* fAidaFile;
    TFile* fBigripsFile;
    TFile* fBrikenFile;

    //! tree to be read
    TTree* ftrAIDA;

    TTree* ftrAIDABeta;
    TTree* ftrAIDAIon;



    TTree* ftrBigrips;
    TTree* ftrNeutron;
    TTree* ftrGamma;
    TTree* ftrAnc;

    //! number of enries in each tree
    Long64_t fnentriesAIDA;
    Long64_t fnentriesAIDABeta;
    Long64_t fnentriesAIDAIon;
    Long64_t fnentriesBigrips;
    Long64_t fnentriesNeutron;
    Long64_t fnentriesGamma;
    Long64_t fnentriesAnc;
    Long64_t fnentriesF11Anc;


    //! data read from stream

    AIDASimpleStruct* faida;

    //AIDA* faidaBeta;
    //AIDA* faidaIon;

    TreeData* fbigrips;
    CloverHit* fclover;
    BELENHit* fneutron;
    BELENHit* fanc;


    //! Declare maps and their iterators (timestamp+entry)
    std::multimap < unsigned long long, unsigned int> fbigripsMap;
    std::multimap < unsigned long long, AIDASimpleStruct* > faidaBetaMap;
    std::multimap < unsigned long long, AIDASimpleStruct* > faidaIonMap;
    std::multimap < unsigned long long, pair<ImplantCorrelationVector*, AIDASimpleStruct* > > faidaImplantMap;
    std::multimap < unsigned long long, BELENHit*> fhe3Map;
    std::multimap < unsigned long long, unsigned int> fcloverMap;
    std::multimap < unsigned long long, unsigned int> fancMap;
    std::multimap < unsigned long long, unsigned int>::iterator fbigripsMap_it;
    std::multimap < unsigned long long, AIDASimpleStruct* >::iterator faidaBetaMap_it;
    std::multimap < unsigned long long, AIDASimpleStruct* >::iterator faidaIonMap_it;
    std::multimap < unsigned long long, pair<ImplantCorrelationVector*, AIDASimpleStruct* > >::iterator faidaImplantMap_it;
    std::multimap < unsigned long long, BELENHit*>::iterator fhe3Map_it;
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


    //! all veto map
    std::multimap < unsigned long long, unsigned int> fvetoMap;
    std::multimap < unsigned long long, unsigned int>::iterator fvetoMap_it;


     //! tree to be filled
     TTree* ftreeImplant;
     TTree* ftreeBeta;
     TTree* ftreeNeutron;
     TTree* ftreeAncillary;

     TTree* ftree;

     //! data struct to be filled (for decay builder)
     IonBetaMult* flocalimp;
     TClonesArray* flocalbeta;
     IonBeta* flocalbetaS;

     BELENHit* flocalneutron;
     IonBetaMult* flocalimp2;

     //! simple data
     datatype decay;

     //! timewindows
     unsigned short fmaxmult;
     double fmaxnpixels;
     long long fIonBetaTWlow;
     long long fIonBetaTWup;
     long long fIonPidTWup;
     long long fIonPidTWlow;
     long long fIonNeutronTWup;
     long long fIonNeutronTWlow;
     long long fIonGammaTWup;
     long long fIonGammaTWlow;
     long long fIonAncTWup;
     long long fIonAncTWlow;


     long long fNeuBetaTWup;
     long long fNeuBetaTWlow;
     long long fNeuBetaoffset;

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

     long long fF11LGammaTWlow;
     long long fF11LGammaTWup;

     double fminneue;
     double fmaxneue;
     double fmineneuvetodown;
     double fmineneuf11;
     Int_t nionbetacorr;

     Long64_t ff11vetodeadtime;
     Long64_t fdownstreamvetodeadtime;
     Long64_t ff11vetototaltime;
     Long64_t fdownstreamvetototaltime;

     Long64_t fvetodeadtime;
     Long64_t fvetototaltime;



     //! stuff for PID separation
     Int_t nri;
     Int_t nbinszet,nbinsaoq;
     Double_t zetrange[2];
     Double_t aoqrange[2];
     Int_t enablepid[MaxNRI];
     Int_t enablepid2[MaxNRI];
     TString nameri[MaxNRI];
     TString latexnametri[MaxNRI];
     Double_t parmsri[MaxNRI][7];
     Double_t halflife[MaxNRI];
     TCutG* cutg[MaxNRI];
     TLatex* pidtag[MaxNRI];

     char* finputMerged;
     TFile* fmergedFile;
     TTree* ftrMerged;
     Long64_t fnentriesMerged;

     TTree* ftreeRI[MaxNRI];
     TTree* ftreeimplantRI[MaxNRI];
     TTree* ftreeimplantAll;
     Long64_t fnentriesImp;




     Double_t dtioncut;
     Int_t sumexyrankcut;
     Int_t lightp_nzcut;
     Double_t neutronecut[2];
     Double_t deltaxcut;
     Double_t deltaycut;


     //! gamma calibration provided by JJ
     Double_t fsep[8];
     Double_t flow_offset[8];
     Double_t flow_gain[8];
     Double_t flow_se[8];
     Double_t fhigh_offset[8];
     Double_t fhigh_gain[8];
     Double_t fhigh_se[8];

     Double_t fcgainold[8];
     Double_t fcoffsetold[8];



     //! temp
     TH1F* fh1;

};

#endif // MERGER_H
