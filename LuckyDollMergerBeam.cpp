#include <fstream>
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
#include "TVectorD.h"

#include "Merger.h"

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

    cout << "Belen Reader" << endl;
    int Verbose = 0;
    char* InputBELEN = NULL;
    char* OutFile = NULL;

    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-bl", "BELEN input list of files", &InputBELEN);
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);

    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      return 1;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      return 2;
    }




    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();

    TTree* treeBeam = new TTree("beam","beam");

    //! Read list of files
    string inputfiles_briken[1000];
    ifstream inf_briken(InputBELEN);

    Int_t nfiles_briken;
    inf_briken>>nfiles_briken;
    TVectorD adruntime(nfiles_briken+1);
    adruntime[0] = 0;
    for (Int_t i=0;i<nfiles_briken;i++){
        adruntime[i+1] = 0;
        inf_briken>>inputfiles_briken[i];
        cout<<inputfiles_briken[i]<<endl;
    }
    for (Int_t i=0;i<nfiles_briken;i++){
        Merger *merge=new Merger();
        merge->SetBrikenFile((char*)inputfiles_briken[i].c_str());
        merge->InitBRIKEN();
        merge->ReadBRIKEN();
        merge->BookTreeBRIKEN(treeBeam);
        ofile->cd();
        merge->DoMergeAnc();
        delete merge;
    }
    treeBeam->Write();
    ofile->Close();
    //! Finish----------------
    double time_end = get_time();
    cout << "\nProgram Run time: " << time_end - time_start << " s." << endl;
    timer.Stop();
    cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

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

    cout << "AIDA event builder" << endl;


}
