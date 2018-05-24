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

using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();


typedef struct {
    unsigned long long T; 	 // Calibrated time
    unsigned long long Tfast;
    double E; 	 // Energy
    double EX;
    double EY;
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int nx, ny, nz;// multiplicity of hits in x and y strips, and dssd planes
    unsigned char ID; 	 // Detector type (BigRips, Aida ion, AIDA beta, BELEN, Clovers)
    //** other stuff pending to define **//
} datatype;


const unsigned char IDion = 4;
const unsigned char IDbeta = 5;
const unsigned char IDcorr = 6;

int main(int argc, char* argv[]){

  //! Program start time
  double time_start = get_time();
  TStopwatch timer;
  timer.Start();
  //! Add signal handler
  signal(SIGINT,signalhandler);

  cout << "AIDA event builder" << endl;
  int Verbose = 0;
  //long long int WindowIon = 5000; //time unit: 10 ns
  //long long int WindowBeta = 2500; //time unit: 10 ns
  //long long int TransientTime = 20000;

  long long int WindowDiscriminator = 0;
  long long int Window = 201;

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

  int SumMultCut=10000;

  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
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
    return 1;
  }
  if(CalibrationFile == NULL){
    cout << "No Calibration table given " << endl;
    CalibrationFile = new char[600];
    strcpy(CalibrationFile,"/sssewqewwq/");
    return 1;
  }
  if(CalibrationFileHE == NULL){
    cout << "No Calibration table for high energy is given " << endl;
    CalibrationFileHE = new char[100];
    strcpy(CalibrationFileHE,"dummy2312039290");
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

  datatype aida;
  TTree* tree=new TTree("aida","aida tree ion and beta)");
  if (FillFlag){
  //! Book tree and histograms
      tree->Branch("aida",&aida,"T/l:Tfast/l:E/D:EX/D:EY/D:x/D:y/D:z/D:nx/I:ny/I:nz/I:ID/b");
  }

  
  //! Read list of files
  string inputfiles[1000];
  ifstream inf(InputAIDA);
  Int_t nfiles;
  inf>>nfiles;


  Int_t implantationrate = 0;
  Int_t betarate = 0;
  Int_t corrscalerrate = 0;

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
      //if (i==0) ecutY[i]=250.;
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

  //!  current time offset
  long long currentTSoffset=0;

  for (Int_t i=0;i<nfiles;i++){
      BuildAIDAEvents* evts=new BuildAIDAEvents;
      evts->SetVerbose(Verbose);
      if (GzFlag!=0) evts->SetGzStream();

      evts->SetMappingFile(MappingFile);
      evts->SetThresholdFile(ThresholdFile);
      evts->SetCalibFile(CalibrationFile);
      evts->SetHECalibFile(CalibrationFileHE);
      evts->SetDiscriminatorTimeWindow(WindowDiscriminator);
      evts->SetTimeWindow(Window);
      evts->SetCorrScalerInStream(true);
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

      double local_time_start = get_time();

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
               (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r "<<flush;
            time_last = time_end;
          }                   
          if (evts->IsBETA()) {//real beta
              //!sort the timestamp
              std::multimap<unsigned long long, int> tsvector;
              std::multimap<unsigned long long, int>::iterator tsvector_it;


              //! remove correlation scaler info for now (2018/03/28)
              /*
              for (int i=0;i<evts->GetAIDACORR()->GetMult();i++){
                  tsvector.insert(std::make_pair(evts->GetAIDACORR()->GetHit(i)->GetTimestamp()*ClockResolution,-1));
              }
              */
              for (int i = 0;i<evts->GetAIDABeta()->GetNClusters();i++){
                  tsvector.insert(std::make_pair(evts->GetAIDABeta()->GetCluster(i)->GetTimestamp() * ClockResolution,i));
              }
              //!fill tree according to the time stamp
              for(tsvector_it = tsvector.begin(); tsvector_it != tsvector.end(); tsvector_it++){
                  Int_t i = tsvector_it->second;
                  if (i!=-1){//real beta
                      Double_t ex = evts->GetAIDABeta()->GetCluster(i)->GetXEnergy();
                      Double_t ey = evts->GetAIDABeta()->GetCluster(i)->GetYEnergy();
                      aida.ID = IDbeta;
                      aida.E = (ex+ey)/2;
                      aida.EX = ex;
                      aida.EY = ey;
                      //!If you need time ordered then we need to modify this
                      aida.T = tsvector_it->first;

                      //! adapt data format for new ion-beta position correlation method
                      uint8_t dx = evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionX()-evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionX();
                      uint8_t dy = evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionY()-evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionY();

                      aida.Tfast = dx + 0x10 * dy;
                      //aida.Tfast = (dx&0xF)|(dy<<4&0xF0);


                      aida.x = (double)(evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionX()+evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionX())/2.;
                      aida.y = (double)(evts->GetAIDABeta()->GetCluster(i)->GetMaxHitPositionY()+evts->GetAIDABeta()->GetCluster(i)->GetMinHitPositionY())/2.;
                      aida.z = evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ();

                      //aida.nx=(int)evts->GetAIDABeta()->GetCluster(i)->GetXMultiplicity();
                      //aida.ny=(int)evts->GetAIDABeta()->GetCluster(i)->GetYMultiplicity();
                      //aida.nz = (int)evts->GetAIDABeta()->GetNClustersZi(int(aida.z));

                      aida.nx = (int)evts->GetAIDABeta()->GetMultX((int)aida.z);
                      aida.ny = (int)evts->GetAIDABeta()->GetMultY((int)aida.z);
                      aida.nz = (int)evts->GetAIDABeta()->GetZHitMult();
                      //aida.nz = (int)evts->GetAIDABeta()->GetNClustersZi((int) aida.z);
                  }else{//corrlation scaler
                      aida.ID=IDcorr;
                      aida.T=tsvector_it->first;
                  }

                  if (evts->GetAIDABeta()->GetCluster(i)->GetYEnergy()>100)
                  if (abs((long long)evts->GetAIDABeta()->GetCluster(i)->GetXTimestamp()*ClockResolution-(long long)evts->GetAIDABeta()->GetCluster(i)->GetYTimestamp()*ClockResolution)<5000)
                  if (evts->GetAIDABeta()->GetCluster(i)->GetSumEXYRank()==0)
                  if (evts->GetAIDABeta()->GetMult()<400) //reject pulser events
                  if (FillFlag) tree->Fill();
              }
          }else if (!evts->IsBETA()){//ion
              //!sort the timestamp
              std::multimap<unsigned long long, int> tsvector;
              std::multimap<unsigned long long, int>::iterator tsvector_it;

              for (int i=0;i<evts->GetAIDACORR()->GetMult();i++){
                  tsvector.insert(std::make_pair(evts->GetAIDACORR()->GetHit(i)->GetTimestamp()*ClockResolution,-1));
              }
              for (int i = 0;i<evts->GetAIDAIon()->GetNClusters();i++){
                  if(evts->GetAIDAIon()->GetCluster(i)->GetHitPositionZ()==evts->GetAIDAIon()->GetMaxZ()) tsvector.insert(std::make_pair(evts->GetAIDAIon()->GetCluster(i)->GetTimestamp() * ClockResolution,i));
              }
              //!fill tree according to the time stamp
              for(tsvector_it = tsvector.begin(); tsvector_it != tsvector.end(); tsvector_it++){
                  Int_t i = tsvector_it->second;
                  if (i!=-1){//real beta
                      Double_t ex = evts->GetAIDAIon()->GetCluster(i)->GetXEnergy();
                      Double_t ey = evts->GetAIDAIon()->GetCluster(i)->GetYEnergy();
                      aida.ID = IDion;
                      aida.E = (ex+ey)/2;
                      aida.EX = ex;
                      aida.EY = ey;
                      //!If you need time ordered then we need to modify this
                      aida.T = tsvector_it->first;

                      uint8_t dx = evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionX()-evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionX();
                      uint8_t dy = evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionY()-evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionY();

                      aida.Tfast = dx + 0x10 * dy;
                      //aida.Tfast = (dx&0xF)|(dy<<4&0xF0);

                      aida.x = (double)(evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionX()+evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionX())/2.;
                      aida.y = (double)(evts->GetAIDAIon()->GetCluster(i)->GetMaxHitPositionY()+evts->GetAIDAIon()->GetCluster(i)->GetMinHitPositionY())/2.;


                      //aida.nx=(int)evts->GetAIDAIon()->GetCluster(i)->GetXMultiplicity();
                      //aida.ny=(int)evts->GetAIDAIon()->GetCluster(i)->GetYMultiplicity();
                      //aida.nz = (int)evts->GetAIDAIon()->GetNClustersZi(int(aida.z));

                      aida.nx = (int)evts->GetAIDAIon()->GetMultX((int)aida.z);
                      aida.ny = (int)evts->GetAIDAIon()->GetMultY((int)aida.z);
                      //aida.nz = (int)evts->GetAIDAIon()->GetClustersMultZ();
                      aida.nz = (int)evts->GetAIDAIon()->GetZHitMult();
                  }else{//corrlation scaler
                      aida.ID=IDcorr;
                      aida.T=tsvector_it->first;
                  }
                  if (FillFlag) tree->Fill();
              }
          }else{//pulser (or something else?))
              for (int i=0;i<evts->GetAIDACORR()->GetMult();i++){
                  aida.ID=IDcorr;
                  aida.T=evts->GetAIDACORR()->GetHit(i)->GetTimestamp()*ClockResolution;
                  if (FillFlag) tree->Fill();
              }
          }
          //alway remember this clear for correlation scaler to avoid memory flow!
          if (evts->GetAIDACORR()->GetMult()>0) evts->GetAIDACORR()->Clear();

          //!Get run time
          if (evts->IsBETA()&&start==0) {
              tstart = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
              start++;
          }else if (!evts->IsBETA()&&start==0){
              tstart = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
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

      cout<<"t"<<tend<<"-"<<tstart<<endl;
      runtime[i+1] = (double)((tend-tstart)*ClockResolution)/(double)1e9;
      runtime[0] += runtime[i+1];
      cout<<"Summary for subrun: "<<inputfiles[i]<<endl;
      cout<<evts->GetCurrentBetaEvent()<<" beta events"<<endl;
      cout<<evts->GetCurrentIonEvent()<<" ion events"<<endl;
      cout<<evts->GetAIDAUnpacker()->GetNCorrScaler()<<" corr events"<<endl;
      implantationrate += evts->GetCurrentIonEvent();
      betarate+=evts->GetCurrentBetaEvent();
      corrscalerrate+=evts->GetAIDAUnpacker()->GetNCorrScaler();
      cout<<ttotal<<" all events (beta+ion)"<<endl;
      cout<<evts->GetCurrentPulserEvent()<<" pulser events"<<endl;
      delete evts;
  }
  if (FillFlag){
      tree->Write();
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
  cout<<"Total number of implantation =  "<<implantationrate<<endl;
  cout<<"Total number of beta =  "<<betarate<<endl;
  cout<<"Total number of corrlation scaler =  "<<corrscalerrate<<endl;

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
