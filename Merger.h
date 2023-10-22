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
#include "TH2.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"

#include "fstream"


#define MaxNRI 1000


const Int_t kMaxGamma = 500;
const Int_t kMaxNeutron = 500;

typedef struct {
    ULong64_t evt;
    ULong64_t ts; 	 //timestamp in ns
    Double_t t,x,y,ex,ey,ion_x,ion_y,ion_ex,ion_ey,zet,aoq,beta,deltaxy;//ts diff in ns
    Short_t z,ion_z,multx,multy,multz,ndecay,isbump;

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
    Double_t neub_E[kMaxNeutron];
    Double_t neub_T[kMaxNeutron];//moderation time in us
    Int_t neub_ch[kMaxNeutron];
    Double_t neub_x[kMaxNeutron];
    Double_t neub_y[kMaxNeutron];

} datatype;

typedef struct {
    ULong64_t evt;
    ULong64_t ts; 	 //timestamp in ns
    Double_t ion_t,ion_x,ion_y,ion_ex,ion_ey,zet,aoq,beta;
    Double_t F11L_T,F11L_E,F11R_T,F11R_E,F7_T,F7_E,veto_T,veto_E,de_T,de_E;
    Short_t ion_z,multx,multy,multz;

    Int_t gc_hit;
    Double_t gc_E[kMaxGamma];
    Double_t gc_T[kMaxGamma];//gamma time in ns
    Double_t gc_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc_ch[kMaxGamma];

    Int_t gc1_hit;
    Double_t gc1_E[kMaxGamma];
    Double_t gc1_T[kMaxGamma];//gamma time in ns
    Double_t gc1_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc1_ch[kMaxGamma];

    Int_t gc2_hit;
    Double_t gc2_E[kMaxGamma];
    Double_t gc2_T[kMaxGamma];//gamma time in ns
    Double_t gc2_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc2_ch[kMaxGamma];

    Int_t ab1_hit;
    Double_t ab1_E[kMaxGamma];
    Double_t ab1_T[kMaxGamma];//gamma time in ns
    Double_t ab1_Tslew[kMaxGamma];//gamma time in ns
    Int_t ab1_ch[kMaxGamma];//first hit channel
    Short_t ab1_mult[kMaxGamma];//multiplicity

    Int_t ab2_hit;
    Double_t ab2_E[kMaxGamma];
    Double_t ab2_T[kMaxGamma];//gamma time in ns
    Double_t ab2_Tslew[kMaxGamma];//gamma time in ns
    Int_t ab2_ch[kMaxGamma];//first hit channel
    Short_t ab2_mult[kMaxGamma];//multiplicity

    Int_t neu_hit;
    Double_t neu_E[kMaxNeutron];
    Double_t neu_T[kMaxNeutron];//moderation time in us
    Int_t neu_ch[kMaxNeutron];

} datatypeisomer;

typedef struct{
    Double_t gc_E;
    Double_t gc_T;//gamma time in ns
    Int_t gc_ch;
} gammahit;

typedef struct{
    Double_t ab_E;
    unsigned long long ab_T;//gamma time in ns
    Double_t ab_Tslew;//gamma time in ns
    Int_t ab_ch;//first hit channel
    Short_t ab_mult[4];//multiplicity
} gammaab;


typedef struct {
    double T; 	 // Calibrated time
    double Tcorr; //correlated time
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int type;
    int type2;
    int evt;
} datatypesimulation;

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
    // normal mode of merging
    void ReadBigrips();  
    void ReadAIDA(unsigned int start = 0, unsigned int stop = 0);
    void ReadBRIKEN(unsigned int startN=0, unsigned int stopN=0,unsigned int startG=0, unsigned int stopG=0,unsigned int startA=0, unsigned int stopA=0);
    void DoMergeSingle();    

    // isomer mode
    void ReadAIDAImpOnly(unsigned int start = 0, unsigned int stop = 0);
    void ReadBRIKENAncGammaOnly(unsigned int startN=0, unsigned int stopN=0,unsigned int startG=0, unsigned int stopG=0,unsigned int startA=0, unsigned int stopA=0);
    void DoMergeIsomer();


    void BookTreeSingle(TTree* tree);
    void BookTreeNeutron(TTree* tree);
    void BookTreeImplant(TTree* tree);

    void SetNeutronOffsetTime(long long offset){fNeuBetaoffset=offset;}
    long long GetNeutronOffsetTime(){return fNeuBetaoffset;}
    void SetSumERankCut(Int_t rankcut){sumexyrankcut = rankcut;}   
    Int_t GetSumERankCut(){return sumexyrankcut;}
    void SetXYTDiffCut(Int_t tcut){xytdiffcut = tcut;}
    Int_t GetXYTDiffCut(){return xytdiffcut;}

    //! for pid separation
    void ReadPID(char* pidfile,Int_t ncutpts=20);
    void SetMergedFile(char* mergedfile){finputMerged = mergedfile;}
    void InitPIDSep();
    void BookPIDSepSimpleTree();
    //! for isomer
    void BookIsomerSimpleTree();


    Int_t GetNri(){return nri;}
    TTree* GetTreeRI(Int_t i){if (i<0) return ftreeallRI; else return ftreeRI[i];}
    TTree* GetTreeImpRI(Int_t i){if (i<0) return ftreeimplantAll;return ftreeimplantRI[i];}

    TCutG* GetCUTRI(Int_t i){return cutg[i];}

    void SetOverlapAreaCorr(){fisoverlapareacorr=true;}

    void SeparatePIDFinal();
    void DoSeparatePIDFinalTree();

    void CopyAddbackData(gammaab* ab_src,gammaab* ab_des);
    void DoAddback();

    TH1F* GetHist1(){return fh1;}
    TH1F* GetHist2(){return fh2;}


    Long64_t GetF11VetoDeadtime(){return ff11vetodeadtime;}
    Long64_t GetDownstreamVetoDeadtime(){return fdownstreamvetodeadtime;}
    Long64_t GetFinalVetoDeadtime(){return fvetodeadtime;}

    Long64_t GetF11VetoTotaltime(){return ff11vetototaltime;}
    Long64_t GetDownstreamVetoTotaltime(){return fdownstreamvetototaltime;}
    Long64_t GetFinalVetoTotaltime(){return fvetototaltime;}


    TH2F* GetH2Deadtime(){return fh2deadtime;}
    TH1F* GetH1Deadtime(){return fh1deadtime;}

    TH2F* GetH2Deadtime2(){return fh2deadtime2;}
    TH1F* GetH1Deadtime2(){return fh1deadtime2;}

    TH2F* GetH2Deadtime3(){return fh2deadtime3;}
    TH1F* GetH1Deadtime3(){return fh1deadtime3;}

    TH2F* GetH2Deadtime4(){return fh2deadtime4;}
    TH1F* GetH1Deadtime4(){return fh1deadtime4;}


    TH2F* GetH2D1(){return fh2d1;}
    TH2F* GetH2D2(){return fh2d2;}


    TH1F* GetPulserHists(Int_t id){return fhpulser[id];}
    TH1F* GetPulserHistsAll(Int_t id){return fhpulserall[id];}


    TH1F* GetH1DeadtimePulserChannel(){return fh1dtpulser;}

    Double_t GetTotalTimePulser(){return ftotaltimepulser;}

    void ResetSimpleData();


    //! reset isomer data
    void ResetIsomerData();

    void BookSimulationTree();
    TTree* GetTreeSimIon(){return ionsimtree;}
    TTree* GetTreeSimBeta(){return betasimtree;}
    TTree* GetTreeSimNeutron(){return neutronsimtree;}

    void DisableImplantNoiseFilterTime(){fflagimpnoiserej = false;}

    void BookDeadTimeTree(TTree* treedeadtime);


    //!stuff for slew correction
    void ReadSlewCorr();
    Bool_t isslewcorr;
    Double_t a[8];
    Double_t b[8];
    Double_t c[8];
    Double_t d[8];

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
    std::multimap < unsigned long long, IonBetaMult* > faidaImplantMapFull;
    std::multimap < unsigned long long, BELENHit*> fhe3Map;
    std::multimap < unsigned long long, unsigned int> fcloverMap;
    std::multimap < unsigned long long, gammaab*> faddbackclover1Map;
    std::multimap < unsigned long long, gammaab*> faddbackclover2Map;


    std::multimap < unsigned long long, unsigned int> fancMap;
    std::multimap < unsigned long long, BELENHit*> fdtpulserMap;
    std::multimap < unsigned long long, unsigned int>::iterator fbigripsMap_it;
    std::multimap < unsigned long long, AIDASimpleStruct* >::iterator faidaBetaMap_it;
    std::multimap < unsigned long long, AIDASimpleStruct* >::iterator faidaIonMap_it;
    std::multimap < unsigned long long, pair<ImplantCorrelationVector*, AIDASimpleStruct* > >::iterator faidaImplantMap_it;
    std::multimap < unsigned long long, IonBetaMult* > ::iterator faidaImplantMapFull_it;
    std::multimap < unsigned long long, BELENHit*>::iterator fhe3Map_it;
    std::multimap < unsigned long long, unsigned int>::iterator fcloverMap_it;
    std::multimap < unsigned long long, gammaab*>::iterator faddbackclover1Map_it;
    std::multimap < unsigned long long, gammaab*>::iterator faddbackclover2Map_it;

    std::multimap < unsigned long long, unsigned int>::iterator fancMap_it;
    std::multimap < unsigned long long, BELENHit*>::iterator fdtpulserMap_it;

    std::multimap < unsigned long long, pair<unsigned int, AIDASimpleStruct*> > fcorrbetaMap;
    std::multimap < unsigned long long, unsigned int> fcorrimpMap;



    //! for Nishimura san
    std::multimap < unsigned long long, unsigned int> ff11ancMap;
    std::multimap < unsigned long long, unsigned int> ff11ancMap_it;

    std::multimap < unsigned long long, unsigned int> fF11LMap;
    std::multimap < unsigned long long, unsigned int> fF11RMap;
    std::multimap < unsigned long long, unsigned int> fF11LRMap;
    std::multimap < unsigned long long, unsigned int> fVetoTopMap;
    std::multimap < unsigned long long, unsigned int> fVetoBotMap;
    std::multimap < unsigned long long, unsigned int> fVetoDownMap;
    std::multimap < unsigned long long, unsigned int>::iterator fF11MapL_it;
    std::multimap < unsigned long long, unsigned int>::iterator fF11MapR_it;

    std::multimap < unsigned long long, unsigned int>::iterator fF11MapLR_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoTopMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoBotMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fVetoDownMap_it;

    std::multimap < unsigned long long, unsigned int> fdETopMap;
    std::multimap < unsigned long long, unsigned int> fdEBotMap;
    std::multimap < unsigned long long, unsigned int>::iterator fdETopMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fdEBotMap_it;

    std::multimap < unsigned long long, pair<unsigned int, AIDASimpleStruct*> >::iterator fcorrbetaMap_it;
    std::multimap < unsigned long long, unsigned int>::iterator fcorrimpMap_it;


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

     TClonesArray* flocalimparray;

     IonBeta* flocalbetaS;

     BELENHit* flocalneutron;

     //! simple decay data
     datatype decay;

     //! simple isomer data
     datatypeisomer isomer;

     //! Merger parameters
     bool fisoverlapareacorr;
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

     //! TW for Isomer
     long long fPIDGammaTWlow;
     long long fPIDGammaTWup;
     long long fPIDVetoTWlow;
     long long fPIDVetoTWup;
     long long fPIDF11TWlow;
     long long fPIDF11TWup;
     long long fPIDdETWlow;
     long long fPIDdETWup;

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

     TTree* ftreeallRI;
     TTree* ftreeRI[MaxNRI];
     TTree* ftreeimplantRI[MaxNRI];
     TTree* ftreeimplantAll;
     Long64_t fnentriesImp;




     Double_t dtioncut;
     Int_t sumexyrankcut;
     Int_t xytdiffcut;
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


     TH1F* fhpulser[141];
     TH1F* fhpulserall[141];


     TH2F* fh2deadtime;
     TH1F* fh1deadtime;


     Long64_t ftsbeginpulser;
     Long64_t ftsendpulser;
     Double_t ftotaltimepulser;
     TH2F* fh2deadtime2;
     TH1F* fh1deadtime2;

     TH2F* fh2deadtime3;
     TH1F* fh1deadtime3;
     TH2F* fh2deadtime4;
     TH1F* fh1deadtime4;

     TH1F* fh1dtpulser;

     //! stuff for rejecting noise events associated with implantation
     Double_t fimpnoisefilter_dxy;
     Short_t fimpnoisefilter_dz;
     Long64_t fimpnoisefilter_dt; //unit ns
     Bool_t fflagimpnoiserej; // unit micro-second



     //! tree of deadtime
     Double_t ftubeno;
     Double_t ftotcnt;
     Double_t ftotcntall;
     Double_t fexpcnt;
     Double_t fdtpulcnt;
     TTree* ftreedeadtime;




     //! temp hist
     TH1F* fh1;
     TH1F* fh2;
     TH2F* fh2d1;
     TH2F* fh2d2;

     //! temp tree
     datatypesimulation ionsim;
     datatypesimulation betasim;
     datatypesimulation neutronsim;
     TTree* ionsimtree;
     TTree* betasimtree;
     TTree* neutronsimtree;


};

#endif // MERGER_H

