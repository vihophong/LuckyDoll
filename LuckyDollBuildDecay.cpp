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

#include "BuildDecay.h"

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

    cout << "Build Decay" << endl;
    int Verbose = 0;
    int Mode = 0;
    char* Input = NULL;

    char* OutFile = NULL;

    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-i", "input list of files", &Input);;
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);
    interface->Add("-m", "Mode selection: 0-build decay 1-beta-neutron for debug 2-ion-neutron for debug (default 0)", &Mode);

    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(Input == NULL){
      cout << "No input list of files given " << endl;
      return 1;
    }

    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      return 2;
    }

    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();

    TTree* treeDecay = new TTree("tree","tree");
    //treeDecay->SetMaxTreeSize(1900000000);

    //! Read list of files
    string inputfiles[1000];
    ifstream inf(Input);

    Int_t nfiles;
    inf>>nfiles;

    TVectorD adruntime(nfiles+1);
    adruntime[0] = 0;
    for (Int_t i=0;i<nfiles;i++){
        adruntime[i+1] = 0;
        inf>>inputfiles[i];
        cout<<"input file; "<<inputfiles[i]<<endl;
    }

    for (Int_t i=0;i<nfiles;i++){
        BuildDecay *decay=new BuildDecay();
        decay->SetInputFile((char*)inputfiles[i].c_str());
        decay->SetMode(Mode);
        decay->Init();
        //! for build decay, one need to book tree later!
        decay->BookTree(treeDecay);
        decay->ShowMeYourFile(ofile);
        decay->ReadImplant();
        decay->ReadBeta();
        decay->ReadNeutron();
        if (Mode==0) {
            decay->DoBuildDecay();
        }
        else if (Mode==1) {
            decay->ReadF11Beam();
            decay->DoBuildDecay3();
        }else if(Mode==2){
            decay->ReadF11Beam();
            decay->DoBuildDecay4();
        }
        ofile->cd();
    }
    treeDecay->Write();
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
