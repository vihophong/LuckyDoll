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

    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-bl", "BELEN input list of files", &InputBELEN);
    interface->Add("-a", "AIDA input list of files", &InputAIDA);
    interface->Add("-br", "Bigrips input list of files", &InputBIGRIPS);
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);

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

    //! Read list of files
    string inputfiles_aida[1000];
    string inputfiles_briken[1000];
    string inputfiles_bigrips[1000];
    ifstream inf_aida(InputAIDA);
    ifstream inf_briken(InputBELEN);
    ifstream inf_bigrips(InputBIGRIPS);

    Int_t nfiles_aida;
    inf_aida>>nfiles_aida;
    Int_t nfiles_briken;
    inf_briken>>nfiles_briken;
    Int_t nfiles_bigrips;
    inf_bigrips>>nfiles_bigrips;

    if (nfiles_aida!=nfiles_briken||nfiles_aida!=nfiles_bigrips){
        cerr<<"error! nfiles briken ne nfiles aida or nfiles_bigrips" <<endl;
        return 3;
    }
    TVectorD adruntime(nfiles_aida+1);
    adruntime[0] = 0;
    for (Int_t i=0;i<nfiles_aida;i++){
        adruntime[i+1] = 0;
        inf_aida>>inputfiles_aida[i];
        inf_briken>>inputfiles_briken[i];
        inf_bigrips>>inputfiles_bigrips[i];
        cout<<inputfiles_aida[i]<<"-"<<inputfiles_briken[i]<<"- "<<inputfiles_bigrips[i]<<endl;
    }

    for (Int_t i=0;i<nfiles_aida;i++){
        Merger *merge=new Merger();
        merge->SetAIDAFile((char*)inputfiles_aida[i].c_str());
        merge->SetBrikenFile((char*)inputfiles_briken[i].c_str());
        merge->SetBigripsFile((char*)inputfiles_bigrips[i].c_str());
        merge->Init();
        merge->ReadAIDA();
        merge->ReadBRIKEN();
        merge->ReadBigrips();
        merge->BookTree(treeImplant,treeBeta,treeNeutron,treeBeam);
        ofile->cd();
        merge->DoMergeImp();
        merge->DoMergeBeta();
        merge->DoMergeNeutron();
        merge->DoMergeAnc();
        delete merge;
    }
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
