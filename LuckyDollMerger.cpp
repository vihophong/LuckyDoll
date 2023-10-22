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

    int XyTDiffCut=5000;


    char* InputMerged = NULL;
    char* InputPID = NULL;

    Int_t CorrMethod=0;   

    Int_t ImplantNoiseRejFlag=0;


    Int_t IsIsomerMode=0;


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

    interface->Add("-isomer", "Isomer mode : 0- disable (default), 1-enable", &IsIsomerMode);

    interface->Add("-bnofs", "Beta Neutron offset time (default=-22000)", &BetaNeutronOffset);

    interface->Add("-rcut", "Sum Energy ranking cut (0- largest energy, a large number(100000)- disabled)", &SumEXYRankcut);

    interface->Add("-tcut", "Time diffrence between X and Y cut (5000 ns by default)", &XyTDiffCut);

    interface->Add("-corrmeth", "Spatial correlation method, 0: old method(default), 1: new method", &CorrMethod);


    interface->Add("-impnoiserej", "Time cut to reject low energy noise events, 1 reject , 0 not reject(default)", &ImplantNoiseRejFlag);


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
    TVectorD deadtimecontainer(7);
    for (Int_t i=0;i<7;i++) deadtimecontainer[i]=0;

    if (!(InputAIDA==NULL)&&InputPID==NULL){
        cout<<"Merging with original merged structure"<<endl;
        TTree* tree = new TTree("tree","tree");
        TTree* treeneutron = new TTree("treeneutron","treeneutron");
        TTree* treeimplant = new TTree("treeimplant","treeimplant");

        Merger *merge=new Merger();
        merge->SetAIDAFile(InputAIDA);
        merge->SetBrikenFile(InputBELEN);
        merge->SetBigripsFile(InputBIGRIPS);
        merge->SetNeutronOffsetTime((long long)BetaNeutronOffset);
        merge->Init();
        merge->BookTreeSingle(tree);
        merge->BookTreeNeutron(treeneutron);
        merge->BookTreeImplant(treeimplant);

        cout<<"Set Neutron Beta offset time = "<<merge->GetNeutronOffsetTime()<<endl;
        merge->ReadAIDA();        
        if (InputBIGRIPS!=NULL) merge->ReadBigrips();
        if (InputBELEN!=NULL) merge->ReadBRIKEN();

        if (CorrMethod!=0) {
            cout<<"New correlation method is chosen!"<<endl;
            merge->SetOverlapAreaCorr();
        }
        merge->DoMergeSingle();

        //! get deadtime,total time
        deadtimecontainer[0]=merge->GetFinalVetoDeadtime();
        deadtimecontainer[1]=merge->GetFinalVetoTotaltime();
        deadtimecontainer[2]=merge->GetF11VetoDeadtime();
        deadtimecontainer[3]=merge->GetF11VetoTotaltime();
        deadtimecontainer[4]=merge->GetDownstreamVetoDeadtime();
        deadtimecontainer[5]=merge->GetDownstreamVetoTotaltime();
        deadtimecontainer[6]=merge->GetTotalTimePulser();
        ofile->cd();
        tree->Write();
        treeneutron->Write();
        treeimplant->Write();

        for (Int_t tubeid=0;tubeid<141;tubeid++){
            merge->GetPulserHists(tubeid)->Write();
            merge->GetPulserHistsAll(tubeid)->Write();
        }
        merge->GetHist1()->Write();
        merge->GetHist2()->Write();
        merge->GetH2Deadtime()->Write();
        merge->GetH1Deadtime()->Write();
        merge->GetH2Deadtime2()->Write();
        merge->GetH1Deadtime2()->Write();

        merge->GetH2Deadtime3()->Write();
        merge->GetH1Deadtime3()->Write();

        merge->GetH2Deadtime4()->Write();
        merge->GetH1Deadtime4()->Write();

        merge->GetH1DeadtimePulserChannel()->Write();
        merge->GetH2D1()->Write();
        merge->GetH2D2()->Write();

        deadtimecontainer.Write("deadtime");
        ofile->Close();
        delete merge;
        //! Finish----------------
    }


    
    /*
    if (!(InputMerged==NULL)){
        cout<<"Separate PID and make simple tree"<<endl;
        Merger *merge=new Merger();
        merge->ReadPID(InputPID);
        merge->SetMergedFile(InputMerged);
        merge->InitPIDSep();
        merge->SetSumERankCut(SumEXYRankcut);
        merge->SetXYTDiffCut(XyTDiffCut);

        if (ImplantNoiseRejFlag==0) {
            merge->DisableImplantNoiseFilterTime();
            cout<<"Disabled implant-related noise rejection"<<endl;
        }

        ofile->cd();
        merge->BookPIDSepSimpleTree();
        cout<<"Ranking cut = "<<merge->GetSumERankCut()<<endl;
        cout<<"TDiff cut = "<<merge->GetXYTDiffCut()<<endl;

        merge->BookSimulationTree();
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

        merge->GetHist1()->Write();
        merge->GetHist2()->Write();

        merge->GetTreeSimIon()->Write();
        merge->GetTreeSimBeta()->Write();
        merge->GetTreeSimNeutron()->Write();

        ofile->Close();
        delete merge;
    }
    */    



    //! all together mode (to reduce file size)
    if ((IsIsomerMode==0)&&(!(InputAIDA==NULL))&&(!(InputPID==NULL))){
        cout<<"Merging with final tree"<<endl;
        //TTree* tree = new TTree("tree","tree");
        TTree* treeneutron = new TTree("treeneutron","treeneutron");
        TTree* treeimplant = new TTree("treeimplant","treeimplant");
        TTree* treedeadtime= new TTree("treedeadtime","treedeadtime");
        Merger *merge=new Merger();
        merge->SetAIDAFile(InputAIDA);
        merge->SetBrikenFile(InputBELEN);
        merge->SetBigripsFile(InputBIGRIPS);
        merge->SetNeutronOffsetTime((long long)BetaNeutronOffset);

        merge->Init();
        //merge->BookTreeSingle(tree);
        merge->BookTreeNeutron(treeneutron);
        merge->BookTreeImplant(treeimplant);
        merge->BookDeadTimeTree(treedeadtime);

        cout<<"Set Neutron Beta offset time = "<<merge->GetNeutronOffsetTime()<<endl;
        merge->ReadAIDA();
        if (InputBIGRIPS!=NULL) merge->ReadBigrips();
        if (InputBELEN!=NULL) merge->ReadBRIKEN();

        if (CorrMethod!=0) {
            cout<<"New correlation method is chosen!"<<endl;
            merge->SetOverlapAreaCorr();
        }
        ofile->cd();


        //! stuff from merger reader
        merge->ReadPID(InputPID);
        merge->SetSumERankCut(SumEXYRankcut);
        merge->SetXYTDiffCut(XyTDiffCut);

        if (ImplantNoiseRejFlag==0) {
            merge->DisableImplantNoiseFilterTime();
            cout<<"Disabled implant-related noise rejection"<<endl;
        }
        merge->BookPIDSepSimpleTree();
        cout<<"Ranking cut = "<<merge->GetSumERankCut()<<endl;
        cout<<"TDiff cut = "<<merge->GetXYTDiffCut()<<endl;
        merge->BookSimulationTree();


        merge->DoMergeSingle();
        //! get deadtime,total time
        deadtimecontainer[0]=merge->GetFinalVetoDeadtime();
        deadtimecontainer[1]=merge->GetFinalVetoTotaltime();
        deadtimecontainer[2]=merge->GetF11VetoDeadtime();
        deadtimecontainer[3]=merge->GetF11VetoTotaltime();
        deadtimecontainer[4]=merge->GetDownstreamVetoDeadtime();
        deadtimecontainer[5]=merge->GetDownstreamVetoTotaltime();
        deadtimecontainer[6]=merge->GetTotalTimePulser();
        ofile->cd();

        //! stuff from merger reader
        for (Int_t i=0;i<merge->GetNri();i++){
            merge->GetCUTRI(i)->Write();
        }
        for (Int_t i=0;i<merge->GetNri();i++){
            merge->GetTreeRI(i)->Write();
            merge->GetTreeImpRI(i)->Write();
        }
        merge->GetTreeRI(-1)->Write();
        merge->GetTreeImpRI(-1)->Write();
        //! end of stuff from merger reader


        treeneutron->Write();
        treeimplant->Write();
        treedeadtime->Write();

        for (Int_t tubeid=0;tubeid<141;tubeid++){
            merge->GetPulserHists(tubeid)->Write();
            merge->GetPulserHistsAll(tubeid)->Write();
        }
        merge->GetHist1()->Write();
        merge->GetHist2()->Write();
        merge->GetH2Deadtime()->Write();
        merge->GetH1Deadtime()->Write();
        merge->GetH2Deadtime2()->Write();
        merge->GetH1Deadtime2()->Write();

        merge->GetH2Deadtime3()->Write();
        merge->GetH1Deadtime3()->Write();

        merge->GetH2Deadtime4()->Write();
        merge->GetH1Deadtime4()->Write();

        merge->GetH1DeadtimePulserChannel()->Write();
        merge->GetH2D1()->Write();
        merge->GetH2D2()->Write();

        deadtimecontainer.Write("deadtime");
        ofile->Close();
        delete merge;
        //! Finish----------------
    }        



    //! Isomer Mode
    if ((IsIsomerMode!=0)&&(!(InputAIDA==NULL))&&(!(InputPID==NULL))){
        cout<<"Merging in ISOMER mode"<<endl;
        //TTree* tree = new TTree("tree","tree");
        Merger *merge=new Merger();
        merge->SetAIDAFile(InputAIDA);
        merge->SetBrikenFile(InputBELEN);
        merge->SetBigripsFile(InputBIGRIPS);
        merge->SetNeutronOffsetTime((long long)BetaNeutronOffset);

        merge->Init();

        cout<<"Set Neutron Beta offset time = "<<merge->GetNeutronOffsetTime()<<endl;
        merge->ReadAIDAImpOnly();
        if (InputBIGRIPS!=NULL) merge->ReadBigrips();
        if (InputBELEN!=NULL) merge->ReadBRIKENAncGammaOnly();
        //if (InputBELEN!=NULL) merge->ReadBRIKEN(0,0,0,0,0,0);

        merge->ReadSlewCorr();
        merge->DoAddback();
        ofile->cd();
        //! stuff from merger reader
        merge->ReadPID(InputPID);
        merge->BookIsomerSimpleTree();

        merge->DoMergeIsomer();
        ofile->cd();

        //! stuff from merger reader
        for (Int_t i=0;i<merge->GetNri();i++){
            merge->GetCUTRI(i)->Write();
        }
        for (Int_t i=0;i<merge->GetNri();i++){
            merge->GetTreeRI(i)->Write();
        }
        merge->GetTreeRI(-1)->Write();

        merge->GetHist1()->Write();
        merge->GetHist2()->Write();

        //! end of stuff from merger reader
        ofile->Close();
        delete merge;
        //! Finish----------------
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
