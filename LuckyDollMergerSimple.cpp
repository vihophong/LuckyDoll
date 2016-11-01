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
#include "BELEN.h"
#include "TVectorD.h"
#include "BelenReader.h"
#include "AIDAUnpacker.h"
#include "BuildAIDAEvents.h"
#include "AIDA.h"

using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();

typedef struct {
    double EX;
    double EY;
    double x,y;// number of pixel for AIDA, or number of tube for BELEN
    unsigned short z;
} aidaSimpleStruct;
typedef struct {
    unsigned short E;
    unsigned short id;
    short x,y;
} he3SimpleStruct;
typedef struct {
    double E;
    short id;// number of pixel for AIDA, or number of tube for BELEN
} cloverSimpleStruct;

int main(int argc, char* argv[]){

    //! Program start time
    double time_start = get_time();
    TStopwatch timer;
    timer.Start();
    //! Add signal handler
    signal(SIGINT,signalhandler);

    cout << "Belen Reader" << endl;
    int Verbose = 0;
    int FillFlag = 1;

    long long int WindowIon = 2500; //time unit: 10 ns
    long long int WindowBeta = 2500; //time unit: 10 ns
    long long int WindowDiscriminator = 0;


    long long int TransientTime = 20000;

    char* InputAIDA = NULL;
    char* CalibrationFile = NULL;
    char* ThresholdFile = NULL;
    char* MappingFile = NULL;
    char* ECutFile = NULL;
    char* He3MappingFile = NULL;



    char* InputBELEN = NULL;

    char* OutFile = NULL;

    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-i", "BELEN input list of files", &InputBELEN);
    interface->Add("-a", "AIDA input list of files", &InputAIDA);
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);
    interface->Add("-mapb", "mapping file (he3 position)", &He3MappingFile);

    interface->Add("-wi", "Ion event building window (default: 2500*10ns)", &WindowIon);
    interface->Add("-wb", "Beta event building window (default: 2500*10ns)", &WindowBeta);
    interface->Add("-wd", "Fast Discriminator Scan window (default: 0 i.e no scan for fast discrimination)", &WindowDiscriminator);

    interface->Add("-map", "mapping file for aida", &MappingFile);
    interface->Add("-cal", "calibration file for aida", &CalibrationFile);
    interface->Add("-thr", "threshold file for aida", &ThresholdFile);

    interface->Add("-f", "fill data or not: 1 fill data 0 no fill (default: fill data)", &FillFlag);
    interface->Add("-tt", "aida transient time (default: 20000*10ns)", &TransientTime);
    interface->Add("-ecut", "specify energy cut file", &ECutFile);


    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      return 1;
    }
    if(He3MappingFile == NULL){
      cout << "No Mapping file given, try to use default: He3_mapping.txt" << endl;
      He3MappingFile = new char[500];
      strcpy(He3MappingFile,"He3_mapping.txt");
      //return 1;
    }
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
      return 1;
    }
    if(CalibrationFile == NULL){
      cout << "No Calibration table given " << endl;
      //return 1;
    }
    if(ECutFile == NULL){
      cout << "No Energy cut file given " << endl;
      //return 2;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      return 2;
    }

    //! Declare maps and their iterators
    std::multimap < unsigned long long, aidaSimpleStruct> aidaBetaMap;
    std::multimap < unsigned long long, aidaSimpleStruct> aidaIonMap;
    std::multimap < unsigned long long, he3SimpleStruct> he3Map;
    std::multimap < unsigned long long, cloverSimpleStruct> cloverMap;
    std::multimap < unsigned long long, aidaSimpleStruct>::iterator aidaBetaMap_it;
    std::multimap < unsigned long long, aidaSimpleStruct>::iterator aidaIonMap_it;
    std::multimap < unsigned long long, he3SimpleStruct>::iterator he3Map_it;
    std::multimap < unsigned long long, cloverSimpleStruct>::iterator cloverMap_it;



    //! READ AIDA
    //! *****************************

    cout<<"Start reading AIDA"<<endl;

    //! Read list of files
    string inputfiles_aida[1000];
    ifstream inf_aida(InputAIDA);
    Int_t nfiles_aida;
    inf_aida>>nfiles_aida;

    Int_t implantationrate = 0;

    TVectorD adruntime(nfiles_aida+1);
    cout<<nfiles_aida<<endl;
    adruntime[0] = 0;
    for (Int_t i=0;i<nfiles_aida;i++){
        adruntime[i+1] = 0;
        inf_aida>>inputfiles_aida[i];
        cout<<inputfiles_aida[i]<<endl;
    }


    Double_t ecutX[6];
    Double_t ecutY[6];
    for(Int_t i=0;i<NumDSSD;i++){
        ecutX[i]=0.;
        ecutY[i]=0.;
        if (i==0) ecutY[i]=250.;
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



    for (Int_t i=0;i<nfiles_aida;i++){
        BuildAIDAEvents* evts=new BuildAIDAEvents;
        evts->SetVerbose(Verbose);
        evts->SetMappingFile(MappingFile);
        evts->SetThresholdFile(ThresholdFile);
        evts->SetCalibFile(CalibrationFile);
        evts->SetDiscriminatorTimeWindow(WindowDiscriminator);
        evts->SetAIDATransientTime(TransientTime);
        evts->SetEventWindowION(WindowIon);
        evts->SetEventWindowBETA(WindowBeta);
        evts->SetPulserInStream(false);
        evts->SetSumEXCut(ecutX);
        evts->SetSumEYCut(ecutY);
        evts->Init((char*)inputfiles_aida[i].c_str());
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
                cout << inputfiles_aida[i] << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" <<
                        (Float_t)ctr/(time_end - local_time_start) << " blocks/s " <<
                        (Float_t)nevtbeta/(time_end - local_time_start) <<" betas/s  "<<
                        (Float_t)nevtion/(time_end - local_time_start) <<" ions/s "<<
                        (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go | ";
                if (evts->IsBETA()) cout<<"beta: x"<<evts->GetAIDABeta()->GetCluster(0)->GetHitPositionX()<<"y"<<
                                          evts->GetAIDABeta()->GetCluster(0)->GetHitPositionY()<<"\r"<< flush;
                else cout<<"ion: x"<<evts->GetAIDAIon()->GetCluster(0)->GetHitPositionX()<<"y"<<
                           evts->GetAIDAIon()->GetCluster(0)->GetHitPositionY()<<"\r"<<flush;
                time_last = time_end;
            }

            if (evts->IsBETA()) {
                //!sort the timestamp
                std::multimap<unsigned long long, int> tsvector;
                std::multimap<unsigned long long, int>::iterator tsvector_it;

                for (int i = 0;i<evts->GetAIDABeta()->GetNClusters();i++){
                    tsvector.insert(std::make_pair(evts->GetAIDABeta()->GetCluster(i)->GetTimestamp() * ClockResolution,i));
                }
                //!fill tree according to the time stamp
                int ii=0;
                for(tsvector_it = tsvector.begin(); tsvector_it != tsvector.end(); tsvector_it++){
                    aidaSimpleStruct aidabeta;
                    Int_t i = tsvector_it->second;
                    Double_t ex = evts->GetAIDABeta()->GetCluster(i)->GetXEnergy();
                    Double_t ey = evts->GetAIDABeta()->GetCluster(i)->GetYEnergy();
                    aidabeta.EX = round(ex);
                    aidabeta.EY = round(ey);
                    //!If you need time ordered then we need to modify this
                    aidabeta.x = round(evts->GetAIDABeta()->GetCluster(i)->GetHitPositionX());
                    aidabeta.y = round(evts->GetAIDABeta()->GetCluster(i)->GetHitPositionY());
                    aidabeta.z = evts->GetAIDABeta()->GetCluster(i)->GetHitPositionZ();
                    aidaBetaMap.insert(std::make_pair(tsvector_it->first,aidabeta));
                    //cout<<tsvector_it->first<<endl;
                    ii++;
                }

            }
            else if (!evts->IsBETA()){
                //! so far only beta!
                int lastclusterID = evts->GetAIDAIon()->GetNClusters()-1;

                if (evts->GetAIDAIon()->GetCluster(lastclusterID)->GetHitPositionZ()==evts->GetAIDAIon()->GetMaxZ()){

                    Double_t ex = evts->GetAIDAIon()->GetCluster(lastclusterID)->GetXEnergy();
                    Double_t ey = evts->GetAIDAIon()->GetCluster(lastclusterID)->GetYEnergy();
                    aidaSimpleStruct aidaion;
                    aidaion.EX = round(ex);
                    aidaion.EY = round(ey);

                    //!If you need time ordered then we need to modify this
                    aidaion.x = round(evts->GetAIDAIon()->GetCluster(lastclusterID)->GetHitPositionX());
                    aidaion.y = round(evts->GetAIDAIon()->GetCluster(lastclusterID)->GetHitPositionY());
                    aidaion.z = evts->GetAIDAIon()->GetCluster(lastclusterID)->GetHitPositionZ();
                    aidaIonMap.insert(std::make_pair(evts->GetAIDAIon()->GetCluster(lastclusterID)->GetTimestamp() * ClockResolution,aidaion));
                }else{
                    cout<<"Somethings wrong with clustering?"<<endl;
                }

            }

            //!Get run time
            if (evts->IsBETA()&&start==0) {
                tstart = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
                start++;
            }else if (!evts->IsBETA()&&start==0){
                tstart = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
                start++;
            }
            if(signal_received){
                break;
            }
        }
        if (evts->IsBETA()) {
            tend = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
        }else if (!evts->IsBETA()){
            tend = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
        }
        adruntime[i+1] = (double)((tend-tstart)*ClockResolution)/(double)1e9;
        adruntime[0] += adruntime[i+1];
        cout<<evts->GetCurrentBetaEvent()<<" beta events"<<endl;
        cout<<evts->GetCurrentIonEvent()<<" ion events"<<endl;
        implantationrate += evts->GetCurrentIonEvent();
        cout<<ttotal<<" all events (beta+ion)"<<endl;
        cout<<evts->GetCurrentPulserEvent()<<" pulser events"<<endl;
        delete evts;
    }






    //! READ BELEN
    //! *****************************
    //! Read list of files
    cout<<"Start reading BELEN"<<endl;
    string inputfiles[1000];
    ifstream inf(InputBELEN);
    Int_t nfiles;
    inf>>nfiles;

    TVectorD runtime(nfiles+1);
    runtime[0] = 0;
    for (Int_t i=0;i<nfiles;i++){
        runtime[i+1] = 0;
        inf>>inputfiles[i];
        cout<<inputfiles[i]<<endl;
    }

    for (Int_t i=0;i<nfiles;i++){
        BelenReader* blrd = new BelenReader;
        blrd->SetMapping(MappingFile);
        //if (FillFlag) blrd->BookTree(treeneutron,treegamma,treeanc);
        blrd->Init((char*)inputfiles[i].c_str());


        int ctr=0;
        int total = blrd->GetNEvents();
        long long tstart = 0;
        long long tend = 0;

        double time_last = (double) get_time();
        int start=0;
        double local_time_start = get_time();

        //! event loop

        while(blrd->GetNextEvent()){
            ctr=blrd->GetCurrentEvent();
            if(ctr%1000000 == 0){
                int nevtneutron = blrd->GetCurrentNeutronEvent();
                int nevtgamma = blrd->GetCurrentGammaEvent();
                int nevtanc = blrd->GetCurrentAncEvent();

                double time_end = get_time();
                //! get time
                cout << inputfiles[i] << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done  " <<
                  (Float_t)(ctr/1000)/(time_end - local_time_start) << "k events/s " <<
                  (Float_t)(nevtneutron/1000)/(time_end - local_time_start) <<"k neutrons/s  "<<
                   (Float_t)(nevtgamma/1000)/(time_end - local_time_start) <<"k gammas/s "<<
                    (Float_t)(nevtanc/1000)/(time_end - local_time_start) <<"k linhtinhs/s "<<
                   (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r "<<flush;
                time_last = time_end;
            }

            //! Get things
            if (blrd->GetType() == 1){//for neutron
                he3SimpleStruct he3simple;
                he3simple.E = blrd->GetEnergy();
                he3simple.id = blrd->GetId();
                he3simple.x = blrd->GetNeutron()->GetPosition().x();
                he3simple.y = blrd->GetNeutron()->GetPosition().y();
                //he3Map.insert(std::make_pair(blrd->GetTimestamp()*20,he3simple));
            }else if (blrd->GetType()==2){//for clover
                cloverSimpleStruct clvsimple;
                clvsimple.E = blrd->GetEnergy();
                clvsimple.id = blrd->GetId();
                cloverMap.insert(std::make_pair(blrd->GetTimestamp()*20,clvsimple));
            }

            //!Get run time
            if (start==0) {
                tstart = blrd->GetTimestamp();
                start++;
            }

            if(signal_received){
                blrd->CloseReader();
              break;
            }
        }

        tend = blrd->GetTimestamp();
        runtime[i+1] = (double) (tend-tstart)/(double)1e9;
        runtime[0] += runtime[i+1];
        cout<<"All Hits= "<<blrd->GetCurrentEvent()<<endl;
        cout<<"  Neutron Hits= "<<blrd->GetCurrentNeutronEvent()<<endl;
        cout<<"  Gamma Hits= "<<blrd->GetCurrentGammaEvent()<<endl;
        cout<<"  Anc Hits= "<<blrd->GetCurrentAncEvent()<<endl;
        blrd->CloseReader();
        delete blrd;

    }

    /*
    //! read old aida files
    //!ion
    std::ifstream inpf("iontable.txt",std::ios::in);
    unsigned long long ts;
    unsigned short z;
    double x,y;
    while (inpf.good()){
        inpf>>ts>>z>>x>>y;
        ts=ts*10;
        aidaSimpleStruct aidaion;
        aidaion.x = x;
        aidaion.y = y;
        aidaion.z = z;
        aidaIonMap.insert(std::make_pair(ts,aidaion));
    }
    //!beta
    std::ifstream inpfb("betatable.txt",std::ios::in);
    while (inpfb.good()){
        inpfb>>ts>>z>>x>>y;
        ts=ts*10;
        aidaSimpleStruct aidabeta;
        aidabeta.x = x;
        aidabeta.y = y;
        aidabeta.z = z;
        aidaBetaMap.insert(std::make_pair(ts,aidabeta));
    }
    */


    cout<<"belen neutron map size = "<<he3Map.size()<<endl;
    cout<<"clover gamma map size = "<<cloverMap.size()<<endl;
    cout<<"aida beta map size = "<<aidaBetaMap.size()<<endl;
    cout<<"aida ion map size = "<<aidaIonMap.size()<<endl;

    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();
    TH1I* tshist = new TH1I("histts","histts",200,-50e3,50e3);
    Long64_t decay_time = 0;
    Long64_t morderate_time = 0;

    unsigned long long time_window1 = 50e3;
    unsigned long long time_window2 = 50e3;
    Long64_t jentry=0;
    Long64_t nionwdecay=0;
    Long64_t nbetawneutron=0;
    Long64_t nbetawgamma=0;
    Long64_t totalsize=aidaBetaMap.size();
    Long64_t check_time = 0;
    Long64_t betats,ionts,he3ts,gammats;



    //! gamma-beta correlation

    for (aidaBetaMap_it=aidaBetaMap.begin();aidaBetaMap_it!=aidaBetaMap.end();aidaBetaMap_it++){
        Double_t percent_complete=(Double_t)jentry/(Double_t)totalsize*100;
        if (jentry%1000==0) cout<<"Finished "<<percent_complete<<" %"<<"| Number of good implantations/total implantations= "<<nbetawgamma<<"/"<<jentry<<endl;
        Long64_t ts1 = (Long64_t)aidaBetaMap_it->first - (Long64_t)time_window1; //100us as dead time
        Long64_t ts2 = (Long64_t)aidaBetaMap_it->first + (Long64_t)time_window2;
        betats = (Long64_t)aidaBetaMap_it->first;
        check_time =  0;
        cloverMap_it = cloverMap.lower_bound(ts1);
        Long64_t ncorr = 0;
        while(cloverMap_it!=cloverMap.end()&&cloverMap_it->first<ts2){
            gammats = (Long64_t)cloverMap_it->first;
            cloverSimpleStruct gamma = cloverMap_it->second;
            if(gammats!=check_time){
                morderate_time  = gammats - betats;
                check_time = gammats;
                tshist->Fill(morderate_time);
                ncorr++;
            }
            cloverMap_it++;
        }
        if (ncorr>0) nbetawgamma++;
        jentry++;
    }


    /*
    //! neutron-beta correlation

    for (aidaBetaMap_it=aidaBetaMap.begin();aidaBetaMap_it!=aidaBetaMap.end();aidaBetaMap_it++){
        Double_t percent_complete=(Double_t)jentry/(Double_t)totalsize*100;
        if (jentry%1000==0) cout<<"Finished "<<percent_complete<<" %"<<"| Number of good implantations/total implantations= "<<nbetawneutron<<"/"<<jentry<<endl;
        Long64_t ts1 = (Long64_t)aidaBetaMap_it->first - (Long64_t)time_window1; //100us as dead time
        Long64_t ts2 = (Long64_t)aidaBetaMap_it->first + (Long64_t)time_window2;
        betats = (Long64_t)aidaBetaMap_it->first;
        check_time =  0;
        he3Map_it = he3Map.lower_bound(ts1);
        Long64_t ncorr = 0;
        while(he3Map_it!=he3Map.end()&&he3Map_it->first<ts2){
            he3ts = (Long64_t)he3Map_it->first;
            he3SimpleStruct he3 = he3Map_it->second;
            if(he3ts!=check_time){
                morderate_time  = he3ts - betats;
                check_time = he3ts;
                tshist->Fill(morderate_time);
                ncorr++;
            }
            he3Map_it++;
        }
        if (ncorr>0) nbetawneutron++;
        jentry++;
    }
    */

    /*
    //! beta-ion correlation
    double deltaxy = 0.;
    for (aidaIonMap_it=aidaIonMap.begin();aidaIonMap_it!=aidaIonMap.end();aidaIonMap_it++){
        Double_t percent_complete=(Double_t)jentry/(Double_t)totalsize*100;
        if (jentry%1000==0) cout<<"Finished "<<percent_complete<<" %"<<"| Number of good implantations/total implantations= "<<nionwdecay<<"/"<<jentry<<endl;

        Long64_t ts1 = (Long64_t)aidaIonMap_it->first - (Long64_t)time_window1; //100us as dead time
        Long64_t ts2 = (Long64_t)aidaIonMap_it->first + (Long64_t)time_window2;
        ionts = (Long64_t)aidaIonMap_it->first;
        aidaSimpleStruct aidaion = aidaIonMap_it->second;

        check_time =  0;
        aidaBetaMap_it = aidaBetaMap.lower_bound(ts1);
        Long64_t ncorr = 0;
        while(aidaBetaMap_it!=aidaBetaMap.end()&&aidaBetaMap_it->first<ts2){
            betats = aidaBetaMap_it->first;
            aidaSimpleStruct aidabeta = aidaBetaMap_it->second;
            if(betats!=check_time&&aidaion.z==aidabeta.z&&aidaion.x>=0&&aidaion.y>=0){
                deltaxy = sqrt((aidaion.x-aidabeta.x)*(aidaion.x-aidabeta.x)+(aidaion.y-aidabeta.y)*(aidaion.y-aidabeta.y));
                if (deltaxy>2.){
                    aidaBetaMap_it++;
                    continue;
                }
                decay_time = (Long64_t)aidaBetaMap_it->first - (Long64_t)aidaIonMap_it->first;
                check_time = betats;
                tshist->Fill(decay_time);

                ncorr++;
            }
            aidaBetaMap_it++;
        }
        if (ncorr>0) nionwdecay++;
        jentry++;
    }
    */

    tshist->Write();


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
