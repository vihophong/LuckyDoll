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

  //! Book tree and histograms
  TTree* treeion=new TTree("ion","tree ion");
  TTree* treebeta=new TTree("beta","tree beta");
  TTree* treepulser=new TTree("pulser","tree pulser");

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
      evts->SetDiscriminatorTimeWindow(WindowDiscriminator);
      evts->SetTimeWindow(Window);
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
               (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r"<<flush;
            time_last = time_end;
          }



          //! special for test
          if (evts->IsBETA()){
              evts->GetBetaTree()->Fill();
              /*
              if (evts->GetAIDABeta()->GetMult()<128){
                  tshist->Fill(evts->GetAIDABeta()->GetTimestamp()-ts_prev_ion_hit);
                  if ((evts->GetAIDABeta()->GetTimestamp()-ts_prev_ion_hit)>4000000){
                      evts->GetAIDABeta()->ClearAllHits();
                      if (evts->GetAIDABeta()->GetNClusters()>0) evts->GetBetaTree()->Fill();
                  }
              }
              */
          }else{
              /*
              if (evts->GetAIDAIon()->GetMult()>0) {
                  //if ((evts->GetAIDAIon()->GetHit(0)->GetTimestamp()-ts_prev_ion_hit)>400000) tdeadtimesum+=1;
                  ts_prev_ion_hit=evts->GetAIDAIon()->GetHit(evts->GetAIDAIon()->GetMult()-1)->GetTimestamp();
              }else{
                  cout<<"sth wrong"<<endl;
              }
              //evts->GetAIDAIon()->ClearAllHits();
              //if (evts->GetAIDAIon()->GetNClusters()>0) evts->GetIonTree()->Fill();
              */
              evts->GetIonTree()->Fill();
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
              //nion++;
              if(evts->GetAIDAIon()->GetHits().size()>0) tend = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
          }
          if(signal_received){
            break;
          }
          //if (nion>70000) goto l1;
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
  }
 l1:
  if (ts_prev_beta_last-ts_prev_ion_hit>100000) tdeadtimesum+=1;
  runtime[0]=tdeadtimesum;
  if (FillFlag){
      //tshist->Write();
      treeion->Write();
      treebeta->Write();
      treepulser->Write();
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
