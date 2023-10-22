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
    char* InputAIDA = NULL;
    char* InputBELEN = NULL;
    char* InputBIGRIPS = NULL;

    char* OutFile = NULL;

    int AIbeg = 0;
    int AIend = 0;
    int ABbeg = 0;
    int ABend = 0;
    int BNbeg = 0;
    int BNend = 0;
    int BGbeg = 0;
    int BGend = 0;
    int BAbeg = 0;
    int BAend = 0;


    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-bl", "BELEN input file", &InputBELEN);
    interface->Add("-a", "AIDA input file", &InputAIDA);
    interface->Add("-br", "Bigrips input file", &InputBIGRIPS);
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);
    interface->Add("-aibeg", "start entry for AIDA ION", &AIbeg);
    interface->Add("-aiend", "stop entry for AIDA ION", &AIend);
    interface->Add("-abbeg", "start entry for AIDA BETA", &ABbeg);
    interface->Add("-abend", "stop entry for AIDA BETA", &ABend);
    interface->Add("-bnbeg", "start entry for neutron in BRIKEN", &BNbeg);
    interface->Add("-bnend", "stop entry for neutron in BRIKEN", &BNend);
    interface->Add("-bgbeg", "start entry for gamma in BRIKEN", &BGbeg);
    interface->Add("-bgend", "stop entry for gamma in BRIKEN", &BGend);
    interface->Add("-babeg", "start entry for anc in BRIKEN", &BAbeg);
    interface->Add("-baend", "stop entry for anc in BRIKEN", &BAend);

    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      return 1;
    }
    if(InputAIDA == NULL){
      cout << "No AIDA input list of files given " << endl;
      return 1;
    }
    if(InputBIGRIPS == NULL){
      cout << "No Bigrips input list of files given " << endl;
      return 1;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      return 2;
    }




    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();

    TTree* treeImplant = new TTree("implant","implant");
    TTree* treeBeta = new TTree("beta","beta");
    TTree* treeNeutron = new TTree("neutron","neutron");
    TTree* treeBeam = new TTree("beam","beam");

    Merger *merge=new Merger();
    merge->SetAIDAFile(InputAIDA);
    merge->SetBrikenFile(InputBELEN);
    merge->SetBigripsFile(InputBIGRIPS);
    merge->Init();
    merge->ReadAIDA(AIbeg,AIend);
    merge->ReadBRIKEN(BNbeg,BNend,BGbeg,BGend,BAbeg,BAend);
    merge->ReadBigrips();
    merge->BookTree(treeImplant,treeBeta,treeNeutron,treeBeam);
    ofile->cd();
    merge->DoMergeImp();
    merge->DoMergeBeta();
    merge->DoMergeNeutron();
    merge->DoMergeAnc();
    delete merge;

    treeImplant->Write();
    treeBeta->Write();
    treeNeutron->Write();
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
