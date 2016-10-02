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
  interface->Add("-a", "AIDA input file", &InputAIDA);
  interface->Add("-o", "output file", &OutFile);
  //interface->Add("-c", "cutfile", &CutFile);
  interface->Add("-v", "verbose level", &Verbose);
  interface->Add("-map", "mapping file (default: FEE_table.txt)", &MappingFile);
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

  //! Program start here
  AIDAUnpacker* aidaunpkg = new AIDAUnpacker;
  aidaunpkg->Init(InputAIDA);

  //!Book a tree
  //aidaunpkg->BookTree();
  aidaunpkg->read_mapping(MappingFile);
  cout<<"maa"<<endl;
  if (ThresholdFile!=NULL) aidaunpkg->read_threshold_table(ThresholdFile);
  aidaunpkg->SetVerbose(Verbose);

  cout<<"Trying to get first Sync...."<<endl;
  aidaunpkg->GetFirstSync();


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
            (Float_t)nhits/(time_end - time_start) << " hits/s (average) - "<< nhits <<" hits "<<
            (total-ctr)*(time_end - time_start)/(Float_t)ctr << "s to go \r" << flush;
          time_last = time_end;
      }
      if(signal_received){
        break;
      }
  }
  cout<<"nhits"<< aidaunpkg->GetHitNumber()<<endl;
  //aidaunpkg->GetTree()->Write();
  ofile->Close();

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
