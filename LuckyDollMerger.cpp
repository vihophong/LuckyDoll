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

    int BetaNeutronOffset=-22000;
    int SumEXYRankcut=0;

    char* InputMerged = NULL;
    char* InputPID = NULL;

    Int_t Mode=0;


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

    interface->Add("-i", "Input Merged file", &InputMerged);

    interface->Add("-pid", "Input PID file", &InputPID);

    interface->Add("-m","Make Simple Tree (0- just separated merged file, 1- make simple saparated tree",&Mode);
    interface->Add("-bnofs", "Beta Neutron offset time (default=-22000)", &BetaNeutronOffset);

    interface->Add("-rcut", "Sum Energy ranking cut (0- largest energy, a large number(100000)- disabled)", &SumEXYRankcut);

    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      //return 1;
    }
    if(InputAIDA == NULL){
      cout << "No AIDA input list of files given " << endl;
      //return 1;
    }
    if(InputBIGRIPS == NULL){
      cout << "No Bigrips input list of files given " << endl;
      //return 1;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      //return 2;
    }    




    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();
    TVectorD deadtimecontainer(6);
    for (Int_t i=0;i<6;i++) deadtimecontainer[i]=0;

    if (!(InputAIDA==NULL)){
        cout<<"Merging"<<endl;
        TTree* tree = new TTree("tree","tree");
        TTree* treeneutron = new TTree("treeneutron","treeneutron");
        TTree* treeimplant = new TTree("treeimplant","treeimplant");

        Merger *merge=new Merger();
        merge->SetAIDAFile(InputAIDA);
        merge->SetBrikenFile(InputBELEN);
        merge->SetBigripsFile(InputBIGRIPS);
        merge->Init();
        merge->SetNeutronOffsetTime((long long)BetaNeutronOffset);
        cout<<"Set Neutron Beta offset time = "<<merge->GetNeutronOffsetTime()<<endl;
        merge->ReadAIDA();        
        if (InputBIGRIPS!=NULL) merge->ReadBigrips();
        if (InputBELEN!=NULL) merge->ReadBRIKEN();

        //merge->BookTreeTClone(tree);
        //merge->DoMergeYOnly();

        //merge->BookTreeTClone(tree);
        //merge->DoMergeTClone();

        merge->BookTreeSingle(tree);
        merge->BookTreeNeutron(treeneutron);
        merge->BookTreeImplant(treeimplant);
        merge->DoMergeSingle();

        /*
        merge->ReadAIDA(AIbeg,AIend);
        merge->ReadBRIKEN(BNbeg,BNend,BGbeg,BGend,BAbeg,BAend);
        merge->ReadBigrips();
        merge->BookTree(treeImplant,treeBeta);
        ofile->cd();
        */

        //! get deadtime,total time
        deadtimecontainer[0]=merge->GetFinalVetoDeadtime();
        deadtimecontainer[1]=merge->GetFinalVetoTotaltime();
        deadtimecontainer[2]=merge->GetF11VetoDeadtime();
        deadtimecontainer[3]=merge->GetF11VetoTotaltime();
        deadtimecontainer[4]=merge->GetDownstreamVetoDeadtime();
        deadtimecontainer[5]=merge->GetDownstreamVetoTotaltime();
        ofile->cd();
        tree->Write();
        treeneutron->Write();
        treeimplant->Write();

        merge->GetHist1()->Write();
        deadtimecontainer.Write("deadtime");
        ofile->Close();
        delete merge;
        //! Finish----------------
    }

    if (!(InputMerged==NULL)){
        if (Mode==0){
            cout<<"Separate PID"<<endl;
            Merger *merge=new Merger();
            merge->ReadPID(InputPID);
            merge->SetMergedFile(InputMerged);
            merge->InitPIDSep();

            deadtimecontainer[0]=merge->GetFinalVetoDeadtime();
            deadtimecontainer[1]=merge->GetFinalVetoTotaltime();
            deadtimecontainer[2]=merge->GetF11VetoDeadtime();
            deadtimecontainer[3]=merge->GetF11VetoTotaltime();
            deadtimecontainer[4]=merge->GetDownstreamVetoDeadtime();
            deadtimecontainer[5]=merge->GetDownstreamVetoTotaltime();
            ofile->cd();
            merge->BookPIDSepTree();
            merge->DoSeparatePID();

            for (Int_t i=0;i<merge->GetNri();i++){
                merge->GetCUTRI(i)->Write();
            }
            for (Int_t i=0;i<merge->GetNri();i++){
                merge->GetTreeRI(i)->Write();
            }
            merge->GetTreeRI(-1)->Write();
            deadtimecontainer.Write("deadtime");
            ofile->Close();
            delete merge;
        }else{
            cout<<"Separate PID and make simple tree"<<endl;
            Merger *merge=new Merger();
            merge->ReadPID(InputPID);
            merge->SetMergedFile(InputMerged);
            merge->InitPIDSep();
            merge->SetSumERankCut(SumEXYRankcut);
            cout<<"Ranking cut = "<<merge->GetSumERankCut()<<endl;

            deadtimecontainer[0]=merge->GetFinalVetoDeadtime();
            deadtimecontainer[1]=merge->GetFinalVetoTotaltime();
            deadtimecontainer[2]=merge->GetF11VetoDeadtime();
            deadtimecontainer[3]=merge->GetF11VetoTotaltime();
            deadtimecontainer[4]=merge->GetDownstreamVetoDeadtime();
            deadtimecontainer[5]=merge->GetDownstreamVetoTotaltime();
            ofile->cd();
            merge->BookPIDSepSimpleTree();
            merge->DoSeparatePIDFinalTree();

            for (Int_t i=0;i<merge->GetNri();i++){
                merge->GetCUTRI(i)->Write();
            }
            for (Int_t i=0;i<merge->GetNri();i++){
                merge->GetTreeRI(i)->Write();
                merge->GetTreeImpRI(i)->Write();
            }
            merge->GetTreeRI(-1)->Write();
            merge->GetTreeImpRI(-1)->Write();
            deadtimecontainer.Write("deadtime");
            ofile->Close();
            delete merge;
        }

    }

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
