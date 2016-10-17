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
#include "AIDAUnpacker.h"
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

  cout << "AIDA Unpacker: Convert raw data into root tree files with no events structure (only hits)" << endl;
  int Verbose = 0;
  char* InputAIDA = NULL;
  char* OutFile = NULL;
  char* ThresholdFile = NULL;
  char* MappingFile = NULL;
  //char* CutFile = NULL;
  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-a", "AIDA list of input files", &InputAIDA);
  interface->Add("-o", "output file", &OutFile);
  //interface->Add("-c", "cutfile", &CutFile);
  interface->Add("-v", "verbose level", &Verbose);
  interface->Add("-map", "mapping file", &MappingFile);
  interface->Add("-thr", "threshold file", &ThresholdFile);


  interface->CheckFlags(argc, argv);
  //Complain about missing mandatory arguments

  if(InputAIDA == NULL){
    cout << "No AIDA input file given " << endl;
    return 1;
  }
  if(MappingFile == NULL){
    cout << "No Mapping file given " << endl;
    return 1;
  }
  if(ThresholdFile == NULL){
    cout << "No Threshold table given " << endl;
    //return 1;
  }
  if(OutFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }

  cout<<"output file: "<<OutFile<< endl;
  TFile* ofile = new TFile(OutFile,"recreate");
  ofile->cd();


  //!Book a tree
  TTree* aidatree = new TTree("tree","raw aida tree");


  //! Read list of files
  string inputfiles[1000];
  ifstream inf(InputAIDA);
  Int_t nfiles;
  inf>>nfiles;

  TVectorD runtime(nfiles+1);
  runtime[0] = 0;
  for (Int_t i=0;i<nfiles;i++){
      runtime[i+1] = 0;
      inf>>inputfiles[i];
      cout<<inputfiles[i]<<endl;
  }

  TH1F* tsCorrHisto;

  for (Int_t i=0;i<nfiles;i++){
      //! Program start here
      AIDAUnpacker* aidaunpkg = new AIDAUnpacker;
      aidaunpkg->Init((char*)inputfiles[i].c_str());
      aidaunpkg->BookTree(aidatree);
      aidaunpkg->read_mapping(MappingFile);
      if (ThresholdFile!=NULL) aidaunpkg->read_threshold_table(ThresholdFile);
      aidaunpkg->SetVerbose(Verbose);

      cout<<"Trying to get first Sync...."<<endl;
      aidaunpkg->GetFirstSync();

      long long tstart;
      long long tend;
      int start=0;

      double local_time_start = get_time();

      cout<<"\n\n****** Start reading AIDA ******\n\n"<<endl;
      double time_last = get_time();
      int ctr = 0;
      int total = aidaunpkg->GetNBlock();
      while (aidaunpkg->GetNextHit()){
          ctr=aidaunpkg->GetCurrentBlock();
          if(ctr%2000 == 0){
              double time_end = get_time();
              long long nhits=aidaunpkg->GetHitNumber();
              cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" <<
                (Float_t)nhits/(time_end - local_time_start) << " hits/s (average) - "<< nhits <<" hits "<<
                (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r" << flush;
              time_last = time_end;
          }
          if (start==0) {
              tstart = aidaunpkg->GetAIDAraw().extTimestamp*tm_stp_scaler_ratio;
          }
          start++;
          tend = aidaunpkg->GetAIDAraw().extTimestamp*tm_stp_scaler_ratio;
          if(signal_received){
            break;
          }
      }

      runtime[i+1] = runtime[i+1] = (double)((tend-tstart)*ClockResolution)/(double)1e9;
      runtime[0] += runtime[i+1];
      cout<<"nhits"<< aidaunpkg->GetHitNumber()<<endl;
      if (i==0) tsCorrHisto = (TH1F*) aidaunpkg->GetTSCorrHist()->Clone();
      else tsCorrHisto->Add(aidaunpkg->GetTSCorrHist());
      delete aidaunpkg;
  }
  aidatree->Write();
  tsCorrHisto->Write();
  runtime.Write("runtime");
  ofile->Close();

  cout<<"\n**********************SUMMARY**********************\n"<<endl;
  cout<<"Total run length = "<<runtime[0]<< " seconds"<<endl;
  cout<<"Sub runs length"<<endl;
  for (Int_t i=0;i<nfiles;i++){
      cout<<inputfiles[i]<<" - "<<runtime[i+1]<< " seconds"<<endl;
  }

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
