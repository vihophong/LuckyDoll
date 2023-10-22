#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "CommandLineInterface.h"
#include "AIDAUnpackerGz.h"
#include "BuildAIDAEventsNew.h"
#include "AIDA.h"
#include "TVectorD.h"
#include "TH2F.h"

typedef struct {
    unsigned long long T; 	 // Calibrated time
    unsigned long long Tfast;
    double E; 	 // Energy
    double EX;
    double EY;
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int nx, ny, nz;// multiplicity of hits in x and y strips, and dssd planes
    int ID; 	 // Detector type (BigRips, Aida ion, AIDA beta, BELEN, Clovers)
    //** other stuff pending to define **//
    int tw; 	 // Calibrated time
    int dz;//dz for dssd1
    int niz;// strip-high-multiplicity flag

} datatype;

const unsigned char IDion = 4;
const unsigned char IDbeta = 5;
const unsigned char IDcorr = 6;


using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();

int main(int argc, char* argv[]){
  //! Program start time
  double time_start = get_time();
  TStopwatch timer;
  timer.Start();
  //! Add signal handler
  signal(SIGINT,signalhandler);

  cout << "AIDA event builder" << endl;
  int Verbose = 0;
  long long int WindowDiscriminator = 0;
  long long int Window = 200;

  int Mode = 0;

  int FillFlag = 1;
  int GzFlag = 0;
  int RankingModeFlag = 1;

  char* InputAIDA = NULL;
  char* OutFile = NULL;
  char* CalibrationFile = NULL;
  char* ThresholdFile = NULL;
  char* MappingFile = NULL;
  char* ECutFile = NULL;
  double ECorr=-1.;
  char* CalibrationFileHE = NULL;

  int MaxDistanceCut=10;
  int MaxDistanceCutIon=3;

  int SumMultCut=10000;

  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-m", "Mode: 0: simple aida tree, 1: advanced aida tree 2: full aidatree", &Mode);
  interface->Add("-a", "AIDA input list of files", &InputAIDA);
  interface->Add("-o", "output file", &OutFile);
  interface->Add("-w", "Time window between successive DAQ readouts (default: 200*10ns)", &Window);
  interface->Add("-wd", "Fast Discriminator Scan window (default: 0 i.e no scan for fast discrimination)", &WindowDiscriminator);
  interface->Add("-v", "verbose level", &Verbose);

  interface->Add("-map", "mapping file", &MappingFile);
  interface->Add("-cal", "calibration file", &CalibrationFile);
  interface->Add("-hecal", "calibration file for high energy", &CalibrationFileHE);
  interface->Add("-thr", "threshold file", &ThresholdFile);

  interface->Add("-f", "fill data or not: 1 fill data 0 no fill (default: fill data)", &FillFlag);
  interface->Add("-ecut", "specify energy cut file", &ECutFile);
  interface->Add("-ecorr", "specify energy corrleration cut", &ECorr);
  interface->Add("-gz", "input data from gz file: 1 enable 0 disable (default: disable)", &GzFlag);

  interface->Add("-rmode", "Switch on(1) off(0) the position determination based on energy correlation ranking (default:on)", &RankingModeFlag);
  interface->Add("-smult", "DSSD multiplicity cut (default 10000)", &SumMultCut);


  interface->CheckFlags(argc, argv);
  //Complain about missing mandatory arguments

  if(InputAIDA == NULL){
    cout << "No AIDA input list of files given " << endl;
    return 1;
  }
  if(MappingFile == NULL){
    cout << "No Mapping file given " << endl;
    return 1;
  }
  if(ThresholdFile == NULL){
    cout << "No Threshold table given " << endl;
    ThresholdFile = new char[600];
    strcpy(ThresholdFile,"/sssewqewwq/");
    //return 1;
  }
  if(CalibrationFile == NULL){
    cout << "No Calibration table given " << endl;
    CalibrationFile = new char[600];
    strcpy(CalibrationFile,"/sssewqewwq/");
    //return 1;
  }
  if(OutFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }
  if(ECutFile == NULL){
    cout << "No Energy cut file given " << endl;
    //return 2;
  }

  cout<<"output file: "<<OutFile<< endl;

  TFile* ofile = new TFile(OutFile,"recreate");
  ofile->cd();
  TH1F* hisx[NumDSSD];
  TH1F* hisy[NumDSSD];
  for (Int_t i=0;i<NumDSSD;i++){
      hisx[i]=new TH1F(Form("hisX%d",i),Form("hisX%d",i),128,0,128);
      hisy[i]=new TH1F(Form("hisY%d",i),Form("hisY%d",i),128,0,128);
  }
  TH2F* h2tdiff=new TH2F("h2tdiff","h2tdiff",20,0,20,5000,-5000,5000);
  TH1F* h1asicsdiff=new TH1F("h1asicsdiff","h1asicsdiff",40,-20,20);
  TH2F* h2asicsdiff=new TH2F("h2asicsdiff","h2asicsdiff",10,0,10,10,0,10);


  AIDASimpleStruct *aida=new AIDASimpleStruct;
  //! Book tree and histograms
  TTree* treeion=new TTree("ion","tree ion");
  TTree* treebeta=new TTree("beta","tree beta");
  TTree* treepulser=new TTree("pulser","tree pulser");
  TTree* treeaida=new TTree("aida","aida simple tree");

  if (FillFlag) treeaida->Branch("aida",&aida);

  //! Read list of files
  string inputfiles[1000];
  ifstream inf(InputAIDA);
  Int_t nfiles;
  inf>>nfiles;

  Int_t implantationrate = 0;

  TVectorD runtime(nfiles+1);
  runtime[0] = 0;
  for (Int_t i=0;i<nfiles;i++){
      runtime[i+1] = 0;
      inf>>inputfiles[i];
      cout<<inputfiles[i]<<endl;
  }

  Double_t ecutX[6];
  Double_t ecutY[6];
  for(Int_t i=0;i<NumDSSD;i++){
      ecutX[i]=-5000.;
      ecutY[i]=-5000.;
  }

  std::ifstream ecut(ECutFile,std::ios::in);
  Int_t dssd = 0;
  if (!ecut.fail()){
      std::cout<<"Reading energy cut from "<<ECutFile<<std::endl;
      while (ecut.good()){
          ecut>>dssd;
          ecut>>ecutX[dssd]>>ecutY[dssd];
      }
      for(Int_t i=0;i<NumDSSD;i++){
          cout<<"dssd No. "<<i<<" X energy cut = "<<ecutX[i]<<"| Y energy cut  = "<<ecutY[i]<<endl;
      }
  }
  //Int_t nion=0;
  Long64_t ts_prev_ion_hit=0;
  Long64_t ts_prev_beta_last=0;
  Long64_t tdeadtimesum=0;


  //!  current time offset
  long long currentTSoffset=0;

  for (Int_t i=0;i<nfiles;i++){
      BuildAIDAEvents* evts=new BuildAIDAEvents;
      evts->SetVerbose(Verbose);
      if (GzFlag!=0) evts->SetGzStream();
      if (FillFlag) evts->BookTree(treeion,treebeta,treepulser);

      evts->SetMappingFile(MappingFile);
      evts->SetThresholdFile(ThresholdFile);
      evts->SetCalibFile(CalibrationFile);
      evts->SetHECalibFile(CalibrationFileHE);
      evts->SetDiscriminatorTimeWindow(WindowDiscriminator);
      evts->SetTimeWindow(Window);
      evts->SetCorrScalerInStream(false);
      evts->SetPulserInStream(true);
      evts->SetSumEXCut(ecutX);
      evts->SetSumEYCut(ecutY);
      cout<<"Ecorr= "<<ECorr<<endl;
      evts->SetEnergyCorrCut(ECorr);
      evts->SetSumMultiplicityCut(SumMultCut);
      if (RankingModeFlag==0) evts->SetNoCorrRankingMode();
      evts->Init((char*)inputfiles[i].c_str());

      //! add last corr ts from previous run, added for briken experiment
      if (i>0) evts->GetAIDAUnpacker()->SetCurrentCorrTSoffset(currentTSoffset);

      double time_last = (double) get_time();

      int ctr=0;
      int total = evts->GetADNblock();
      int ttotal=0;

      long long tstart;
      long long tend;
      int start=0;

      int ncnt1=0;
      int ncnt2=0;
      int ncnt3=0;


      double local_time_start = get_time();

      int nbeta=0;
      int nion=0;
      //!event loop
      while(evts->GetNextEvent()){
          ttotal++;
          ctr=evts->GetCurrentADBlock();
          if(ctr%1000 == 0){
            int nevtbeta = evts->GetCurrentBetaEvent();
            int nevtion = evts->GetCurrentIonEvent();
            double time_end = get_time();
            cout << inputfiles[i] << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" <<
              (Float_t)ctr/(time_end - local_time_start) << " blocks/s " <<
              (Float_t)nevtbeta/(time_end - local_time_start) <<" betas/s  "<<
               (Float_t)nevtion/(time_end - local_time_start) <<" ions/s "<<
               (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r"<<flush;
            time_last = time_end;
          }
          //! special for test
          if (evts->IsBETA()){
              //!sort the timestamp
              std::multimap<unsigned long long, int> tsvector;
              std::multimap<unsigned long long, int>::iterator tsvector_it;
              tsvector.clear();
              //! cluster XY map
              for (int i = 0;i<evts->GetAIDABeta()->GetNClusters();i++){
                  tsvector.insert(std::make_pair(evts->GetAIDABeta()->GetCluster(i)->GetTimestamp() * ClockResolution,i));
              }



              //! for vetoing purpose
              Double_t maxeveto=-9999.;
              if (evts->GetAIDABeta()->GetMultY(NumDSSD-1)>0){
                  for (int i = 0;i<evts->GetAIDABeta()->GetMult();i++){
                      if (evts->GetAIDABeta()->GetHit(i)->GetZ()==5&&evts->GetAIDABeta()->GetHit(i)->GetXY()>=128){//select Y strips only
                          if (evts->GetAIDABeta()->GetHit(i)->GetEnergy()>maxeveto) maxeveto = evts->GetAIDABeta()->GetHit(i)->GetEnergy();
                      }
                  }
              }
              Int_t zmult=0; for (int i=0;i<NumDSSD;i++) if (evts->GetAIDABeta()->GetMultY(i)>0) zmult++;

              if (Mode==0){
                  //!fill tree according to the time stamp
                  for(tsvector_it = tsvector.begin(); tsvector_it != tsvector.end(); tsvector_it++){
                      Int_t i = tsvector_it->second;
                      Double_t ex = evts->GetAIDABeta()->GetCluster(i)->GetXEnergy();
                      Double_t ey = evts->GetAIDABeta()->GetCluster(i)->GetYEnergy();
                      aida->Clear();
                      aida->SetID(IDbeta);
                      aida->SetEventNumber(nbeta);
                      aida->SetXEnergy(ex);
                      aida->SetYEnergy(ey);

                      aida->SetXTimestamp(evts->GetAIDABeta()->GetCluster(i)->GetXTimestamp()*ClockResolution);
                      aida->SetYTimestamp(evts->GetAIDABeta()->GetCluster(i)->GetYTimestamp()*ClockResolution);
                      aida->SetXMaxTimestamp(evts->GetAIDABeta()->GetCluster(i)->GetXMaxTimestamp()*ClockResolution);
                      aida->SetYMaxTimestamp(evts->GetAIDABeta()->GetCluster(i)->GetYMaxTimestamp()*ClockResolution);

                      aida->SetXMinPos(evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionX());
                      aida->SetXMaxPos(evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionX());
                      aida->SetYMinPos(evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionY());
                      aida->SetYMaxPos(evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionY());

                      aida->SetTimeDifferenceMin(evts->GetAIDABeta()->GetCluster(i)->GetTimeDifferenceMin()*ClockResolution);
                      aida->SetTimeDifferenceMax(evts->GetAIDABeta()->GetCluster(i)->GetTimeDifferenceMax()*ClockResolution);

                      aida->SetTimestamp(tsvector_it->first);
                      //!If you need time ordered then we need to modify this
                      aida->SetHitPosition(evts->GetAIDABeta()->GetCluster(i)->GetHitPositionX(),evts->GetAIDABeta()->GetCluster(i)->GetHitPositionY(),evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ());
                      aida->SetMult(evts->GetAIDABeta()->GetMult());
                      aida->SetXClusterMult(evts->GetAIDABeta()->GetNXClustersZi((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetYClusterMult(evts->GetAIDABeta()->GetNYClustersZi((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetXMult((unsigned short)evts->GetAIDABeta()->GetMultX((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetYMult((unsigned short)evts->GetAIDABeta()->GetMultY((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetZMult(zmult);
                      //aida->SetZMult((unsigned short)evts->GetAIDABeta()->GetZHitMult());
                      //aida->SetZMult((unsigned short)evts->GetAIDABeta()->GetNClustersZi((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetStripMultFlag((unsigned short)evts->GetAIDABeta()->GetStripMultFlag((int)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetTimeWidth((double)(evts->GetAIDABeta()->GetTimeWindow()*ClockResolution)/1000.);
                      aida->SetDZ((unsigned short) evts->GetAIDABeta()->GetDeltaMaxZ());
                      aida->SetRankingFlag(evts->GetAIDABeta()->GetCluster(i)->GetRankingFlag());
                      aida->SetEDiffRankingFlag(evts->GetAIDABeta()->GetCluster(i)->GetEDiffRankingFlag());
                      aida->SetDtIon((double)(evts->GetAIDABeta()->GetCluster(i)->GetTimestamp()-evts->GetAIDABeta()->GetPrevIonTimestamp())*ClockResolution/1000.);
                      aida->SetMinimumDistanceX(evts->GetAIDABeta()->GetMinDistanceX((unsigned short)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetMinimumDistanceY(evts->GetAIDABeta()->GetMinDistanceY((unsigned short)evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetSumEXYRank(evts->GetAIDABeta()->GetCluster(i)->GetSumEXYRank());
                      aida->SetMaxELastDSSD(maxeveto);
                      //if (abs(aida->GetXEnergy()-aida->GetYEnergy())<300)

                      if (aida->GetYEnergy()>100)
                      if (abs((long long)aida->GetXMinTimestamp()-(long long)aida->GetYMinTimestamp())<5000)
                      if (aida->GetSumEXYRank()==0)
                      if (aida->GetMultiplicity()<400) //reject pulser events
                      if (FillFlag) treeaida->Fill();
                  }
                  if (tsvector.size()>0) nbeta++;
              }

              if (Mode==1) evts->GetAIDABeta()->ClearAllHits();
              if (Mode==1&&evts->GetAIDABeta()->GetNClusters()>0) evts->GetBetaTree()->Fill();
              if (Mode==2&&evts->GetAIDABeta()->GetNClusters()>0) evts->GetBetaTree()->Fill();
          }else{
              //!sort the timestamp
              std::multimap<unsigned long long, int> tsvector;
              std::multimap<unsigned long long, int>::iterator tsvector_it;
              tsvector.clear();

              //! cluster XY map
              std::map < Double_t ,Double_t> fxmap;
              std::map < Double_t ,Double_t> fymap;
              std::map < Double_t ,Double_t> ::iterator fxmap_it;
              std::map < Double_t ,Double_t> ::iterator fymap_it;
              for (int i = 0;i<evts->GetAIDAIon()->GetNClusters();i++){
                  if (evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()==evts->GetAIDAIon()->GetMaxZ()) {
                      tsvector.insert(std::make_pair(evts->GetAIDAIon()->GetCluster(i)->GetTimestamp() * ClockResolution,i));
                      Double_t xpos=evts->GetAIDAIon()->GetCluster(i)->GetHitPositionX();
                      Double_t ypos=evts->GetAIDAIon()->GetCluster(i)->GetHitPositionY();
                      fxmap.insert(std::make_pair(xpos,-1.));
                      fymap.insert(std::make_pair(ypos,-1.));
                  }
              }
              if (evts->GetAIDAIon()->GetMaxZ()<5) {
                  ncnt2++;
                  if (evts->GetAIDAIon()->GetMultX((short)evts->GetAIDAIon()->GetMaxZ()+1)>0||evts->GetAIDAIon()->GetMultY((short)evts->GetAIDAIon()->GetMaxZ()+1)>0){
                      ncnt3++;
                  }
              }
              double mindx=9999;
              double mindy=9999;
              double xprev=-1;
              double yprev=-1;
              for (fxmap_it=fxmap.begin();fxmap_it!=fxmap.end();fxmap_it++){
                  if (xprev>0&&((fxmap_it->first-xprev)<mindx)) mindx=fxmap_it->first-xprev;
                  xprev=fxmap_it->first;
              }
              for (fymap_it=fymap.begin();fymap_it!=fymap.end();fymap_it++){
                  if (yprev>0&&((fymap_it->first-yprev)<mindy)) mindy=fymap_it->first-yprev;
                  yprev=fymap_it->first;
              }

              Int_t zmult=0; for (int i=0;i<NumDSSD;i++) if (evts->GetAIDAIon()->GetMultY(i)>0) zmult++;
              if (Mode==0){
                  //!fill tree according to the time stamp
                  for(tsvector_it = tsvector.begin(); tsvector_it != tsvector.end(); tsvector_it++){
                      Int_t i = tsvector_it->second;
                      Double_t ex = evts->GetAIDAIon()->GetCluster(i)->GetXEnergy();
                      Double_t ey = evts->GetAIDAIon()->GetCluster(i)->GetYEnergy();
                      aida->Clear();
                      aida->SetEventNumber(nion);
                      aida->SetID(IDion);
                      aida->SetXEnergy(ex);
                      aida->SetYEnergy(ey);
                      aida->SetTimestamp(tsvector_it->first);


                      aida->SetXTimestamp(evts->GetAIDAIon()->GetCluster(i)->GetXTimestamp()*ClockResolution);
                      aida->SetYTimestamp(evts->GetAIDAIon()->GetCluster(i)->GetYTimestamp()*ClockResolution);
                      aida->SetXMaxTimestamp(evts->GetAIDAIon()->GetCluster(i)->GetXMaxTimestamp()*ClockResolution);
                      aida->SetYMaxTimestamp(evts->GetAIDAIon()->GetCluster(i)->GetYMaxTimestamp()*ClockResolution);

                      aida->SetXMinPos(evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionX());
                      aida->SetXMaxPos(evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionX());
                      aida->SetYMinPos(evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionY());
                      aida->SetYMaxPos(evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionY());

                      aida->SetTimeDifferenceMin(evts->GetAIDAIon()->GetCluster(i)->GetTimeDifferenceMin()*ClockResolution);
                      aida->SetTimeDifferenceMax(evts->GetAIDAIon()->GetCluster(i)->GetTimeDifferenceMax()*ClockResolution);


                      //!If you need time ordered then we need to modify this
                      aida->SetHitPosition(evts->GetAIDAIon()->GetCluster(i)->GetHitPositionX(),evts->GetAIDAIon()->GetCluster(i)->GetHitPositionY(),evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ());

                      aida->SetMult(evts->GetAIDAIon()->GetMult());
                      aida->SetXClusterMult(evts->GetAIDAIon()->GetNXClustersZi((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetYClusterMult(evts->GetAIDAIon()->GetNYClustersZi((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));

                      aida->SetXMult((unsigned short)evts->GetAIDAIon()->GetMultX((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetYMult((unsigned short)evts->GetAIDAIon()->GetMultY((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetZMult(zmult);
                      //aida->SetZMult((unsigned short)evts->GetAIDAIon()->GetZHitMult());
                      //aida->SetZMult((unsigned short)evts->GetAIDAIon()->GetNClustersZi((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetStripMultFlag((unsigned short)evts->GetAIDAIon()->GetStripMultFlag((int)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetTimeWidth((double)(evts->GetAIDAIon()->GetTimeWindow()*ClockResolution)/1000.);
                      aida->SetDZ((unsigned short) evts->GetAIDAIon()->GetDeltaMaxZ());
                      aida->SetRankingFlag(evts->GetAIDAIon()->GetCluster(i)->GetRankingFlag());
                      aida->SetEDiffRankingFlag(evts->GetAIDAIon()->GetCluster(i)->GetEDiffRankingFlag());
                      aida->SetDtIon((double)(evts->GetAIDAIon()->GetCluster(i)->GetTimestamp()-evts->GetAIDAIon()->GetPrevIonTimestamp())*ClockResolution/1000.);
                      aida->SetMinimumDistanceX(evts->GetAIDABeta()->GetMinDistanceX((unsigned short)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetMinimumDistanceY(evts->GetAIDABeta()->GetMinDistanceY((unsigned short)evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()));
                      aida->SetSumEXYRank(evts->GetAIDAIon()->GetCluster(i)->GetSumEXYRank());
                      //if (abs((long long)aida->GetXMinTimestamp()-(long long)aida->GetYMinTimestamp())<5000)
                      //if (abs(aida->GetXEnergy()-aida->GetYEnergy())<300)
                      if (FillFlag) treeaida->Fill();
                  }
                  if (tsvector.size()>0) nion++;

              }
              if (Mode==1) evts->GetAIDAIon()->ClearAllHits();
              if (Mode==1&&evts->GetAIDAIon()->GetNClusters()>0) evts->GetIonTree()->Fill();
              if (Mode==2&&evts->GetAIDAIon()->GetNClusters()>0) evts->GetIonTree()->Fill();
          }
          if (evts->IsBETA()&&Mode==2){
              if (evts->GetAIDABeta()->GetMult()>64) evts->GetPulserTree()->Fill();
          }

          //!Get run time
          if (evts->IsBETA()&&start==0) {
              if(evts->GetAIDABeta()->GetHits().size()>0) tstart = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
              start++;
          }else if (!evts->IsBETA()&&start==0){
              if(evts->GetAIDAIon()->GetHits().size()>0) tstart = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
              start++;
          }
          if (evts->IsBETA()) {
              if(evts->GetAIDABeta()->GetHits().size()>0) tend = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
          }else if (!evts->IsBETA()){
              if(evts->GetAIDAIon()->GetHits().size()>0) tend = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
          }
          if(signal_received){
            break;
          }

      }
      //! Get last corr ts, added for briken experiment
      currentTSoffset=evts->GetAIDAUnpacker()->GetCurrentCorrTSoffset();

      runtime[i+1] = (double)((tend-tstart)*ClockResolution)/(double)1e6;
      runtime[0] += runtime[i+1];
      cout<<evts->GetCurrentBetaEvent()<<" beta events"<<endl;
      cout<<evts->GetCurrentIonEvent()<<" ion events"<<endl;
      implantationrate += evts->GetCurrentIonEvent();
      cout<<ttotal<<" all events (beta+ion)"<<endl;
      cout<<evts->GetCurrentPulserEvent()<<" pulser events"<<endl;
      delete evts;
      cout<<ncnt1<<endl;
      cout<<ncnt2<<endl;
      cout<<ncnt3<<endl;
  }
 l1:
  if (ts_prev_beta_last-ts_prev_ion_hit>100000) tdeadtimesum+=1;
  runtime[0]=tdeadtimesum;

  if (FillFlag){
      //tshist->Write();
      h2tdiff->Write();
      h1asicsdiff->Write();
      h2asicsdiff->Write();
      treeion->Write();
      treebeta->Write();
      treepulser->Write();
      treeaida->Write();
      runtime.Write("runtime");
  }

  ofile->Close();


  cout<<"\n**********************SUMMARY**********************\n"<<endl;
  cout<<"Total run length = "<<runtime[0]<< " seconds"<<endl;
  cout<<"Sub runs length"<<endl;
  for (Int_t i=0;i<nfiles;i++){
      cout<<inputfiles[i]<<" - "<<runtime[i+1]<< " seconds"<<endl;
  }
  cout<<"\nImplatation rate =  "<<(double)implantationrate/runtime[0]<< " cps"<<endl;
  //runtime.Print();
  //! Finish----------------
  double time_end = get_time();
  cout << "\nProgram Run time: " << time_end - time_start << " s." << endl;
  timer.Stop();
  cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
  return 0;
}

void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}
