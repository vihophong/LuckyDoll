#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "CommandLineInterface.h"
#include "AIDAUnpacker.h"
#include "BuildAIDAEventsNew.h"
#include "AIDA.h"
#include "TVectorD.h"

using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();


typedef struct {
    double EX;
    double EY;
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    //** other stuff pending to define **//
} datatype;


const unsigned char IDion = 4;
const unsigned char IDbeta = 5;

int main(int argc, char* argv[]){

  //! Program start time
  double time_start = get_time();
  TStopwatch timer;
  timer.Start();
  //! Add signal handler
  signal(SIGINT,signalhandler);

  cout << "Generate Pulser or 207Bi calibration spectra" << endl;
  int Verbose = 0;
  int GzFlag = 0;
  int Mode = NULL;
  int Mult = 1;
  long long int WindowIon = 2500; //time unit: 10 ns
  long long int WindowBeta = 2500; //time unit: 10 ns
  long long int WindowDiscriminator = 0;

  double nbins = 4000;
  double min = 0;
  double max = 24000;

  long long int TransientTime = 20000;

  char* InputAIDA = NULL;
  char* OutFile = NULL;
  char* CalibrationFile = NULL;
  char* ThresholdFile = NULL;
  char* MappingFile = NULL;
  char* CalibrationFileHE = NULL;
  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-a", "AIDA input list of files", &InputAIDA);
  interface->Add("-o", "output file", &OutFile);
  //interface->Add("-wi", "Ion event building window (default: 2500*10ns)", &WindowIon);
  //interface->Add("-wb", "Beta event building window (default: 2500*10ns)", &WindowBeta);
  interface->Add("-wd", "Fast Discriminator Scan window (default: 0 i.e no scan for fast discrimination)", &WindowDiscriminator);
  interface->Add("-v", "verbose level", &Verbose);

  interface->Add("-map", "mapping file ", &MappingFile);
  interface->Add("-cal", "calibration file", &CalibrationFile);
  interface->Add("-hecal", "calibration file for high energy", &CalibrationFileHE);
  interface->Add("-thr", "threshold file", &ThresholdFile);

   //interface->Add("-tt", "aida transient time (default: 20000*10ns)", &TransientTime);
   interface->Add("-m", " mode selection: 1: pulser spectra 2: source spectra", &Mode);
   interface->Add("-mult", " Multiplicity condition (<=) in X and Y strips (for source spectra mode)", &Mult);

   interface->Add("-nbins","number of histogram bin (default 4000)",&nbins);
   interface->Add("-min","minimum range of X (default 0)",&min);
   interface->Add("-max","maximum range of X (default 24000)",&max);
   interface->Add("-gz", "input data from gz file: 1 enable 0 disable (default: disable)", &GzFlag);

  interface->CheckFlags(argc, argv);
  //Complain about missing mandatory arguments

  if(Mode == NULL){
    cout << "No Mode selection given " << endl;
    return 1;
  }
  if(Mode >4 || Mode <1){
    cout << "Invalid Mode selection " << endl;
    return 1;
  }


  if(Mult == 1){
    cout << "Use default multiplicity < 1 " << endl;
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
    ThresholdFile = new char[100];
    strcpy(ThresholdFile,"dummy2312039290");
    //return 1;
  }
  if(CalibrationFile == NULL){
    cout << "No Calibration table is given " << endl;
    CalibrationFile = new char[100];
    strcpy(CalibrationFile,"dummy2312039290");
    //return 1;
  }
  if(CalibrationFileHE == NULL){
    cout << "No Calibration table for high energy is given " << endl;
    CalibrationFileHE = new char[100];
    strcpy(CalibrationFileHE,"dummy2312039290");
    //return 1;
  }
  if(OutFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }


  cout<<"output file: "<<OutFile<< endl;

  TFile* ofile = new TFile(OutFile,"recreate");
  ofile->cd();

  TTree* tree;
  datatype aida;
  aida.EY=0;
  aida.EX=0;
  aida.x=0;
  aida.y=0;
  aida.z=0;
  if (Mode==3||Mode==4) {
      tree=new TTree("tree","tree");
      tree->Branch("aida",&aida,"EX/D:EY/D:x/D:y/D:z/D");
  }

  //! Book tree and histograms
  TH2F* spectra=new TH2F ("spectra","spectra",NumDSSD*(NumStrX+NumStrY),0,NumDSSD*(NumStrX+NumStrY),nbins,min,max);
  TH2F* calibspectra=new TH2F ("calibspectra","calibspectra",NumDSSD*(NumStrX+NumStrY),0,NumDSSD*(NumStrX+NumStrY),nbins,min,max);


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

  Double_t ecutX[6];
  Double_t ecutY[6];
  for(Int_t i=0;i<NumDSSD;i++){
      ecutX[i]=-5000.;
      ecutY[i]=-5000.;
  }

  for (Int_t i=0;i<nfiles;i++){
      BuildAIDAEvents* evts=new BuildAIDAEvents;
      evts->SetVerbose(Verbose);
      if (GzFlag!=0) evts->SetGzStream();
      evts->SetMappingFile(MappingFile);
      evts->SetThresholdFile(ThresholdFile);
      evts->SetCalibFile(CalibrationFile);
      evts->SetHECalibFile(CalibrationFileHE);
      evts->SetDiscriminatorTimeWindow(WindowDiscriminator);
      //evts->SetAIDATransientTime(TransientTime);
      //evts->SetEventWindowION(WindowIon);
      //evts->SetEventWindowBETA(WindowBeta);
      evts->SetSumEXCut(ecutX);
      evts->SetSumEYCut(ecutY);
      if (Mode == 1) evts->SetPulserInStream(true);
      else evts->SetPulserInStream(false);
      evts->Init((char*)inputfiles[i].c_str());
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
            int nevtpulser = evts->GetCurrentPulserEvent();
            double time_end = get_time();
            cout << inputfiles[i] << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done\t" <<
              (Float_t)ctr/(time_end - local_time_start) << " blocks/s " <<
              (Float_t)nevtbeta/(time_end - local_time_start) <<" betas/s  "<<
                (Float_t)nevtpulser/(time_end - local_time_start) <<" pulsers/s  "<<
               (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r "<<flush;
            time_last = time_end;
          }
          if ((Mode==1) && evts->IsBETA() && evts->GetAIDABeta()->GetMult()>=64) { //pulser mode
              for (Int_t j = 0;j<evts->GetAIDABeta()->GetMult();j++){
                  int adc = evts->GetAIDABeta()->GetHit(j)->GetADC();
                  double energy = evts->GetAIDABeta()->GetHit(j)->GetEnergy();
                  short z =  evts->GetAIDABeta()->GetHit(j)->GetZ();
                  short xy = evts->GetAIDABeta()->GetHit(j)->GetXY();
                  if (xy<NumStrX) {
                      spectra->Fill(xy + z*NumStrX,adc);
                      calibspectra->Fill(xy + z*NumStrX,energy);
                  }else{
                      spectra->Fill((xy - NumStrX) + z*NumStrY + NumDSSD*NumStrX,adc);
                      calibspectra->Fill((xy - NumStrX) + z*NumStrY + NumDSSD*NumStrX,energy);
                  }
              }
          }
          if ((Mode==2) && evts->IsBETA() && evts->GetAIDABeta()->GetMult()<64) { //beta mode
              for (Int_t j = 0;j<evts->GetAIDABeta()->GetMult();j++){
                  int adc = evts->GetAIDABeta()->GetHit(j)->GetADC();
                  double energy = evts->GetAIDABeta()->GetHit(j)->GetEnergy();
                  short z =  evts->GetAIDABeta()->GetHit(j)->GetZ();
                  short xy = evts->GetAIDABeta()->GetHit(j)->GetXY();
                  if (xy<NumStrX) {
                      if (evts->GetAIDABeta()->GetMultX(z) <= (unsigned short) Mult){
                          spectra->Fill(xy + z*NumStrX,adc);
                          calibspectra->Fill(xy + z*NumStrX,energy);
                      }
                  }else{
                      if (evts->GetAIDABeta()->GetMultY(z) <= (unsigned short) Mult){
                          spectra->Fill((xy - NumStrX) + z*NumStrY + NumDSSD*NumStrX,adc);
                          calibspectra->Fill((xy - NumStrX) + z*NumStrY + NumDSSD*NumStrX,energy);
                      }
                  }
              }
          }
          if ((Mode==3) && !evts->IsBETA()){//generate tree files for implantation calibration
              for (unsigned short jj=0;jj<evts->GetAIDAIon()->GetClusters().size();jj++){
                  if (evts->GetAIDAIon()->GetCluster(jj)->GetXMultiplicity()==1&&evts->GetAIDAIon()->GetCluster(jj)->GetYMultiplicity()==1){
                      aida.z=evts->GetAIDAIon()->GetCluster(jj)->GetHitPositionZ();
                      aida.y=evts->GetAIDAIon()->GetCluster(jj)->GetHitPositionY();
                      aida.x=evts->GetAIDAIon()->GetCluster(jj)->GetHitPositionX();
                      aida.EX=evts->GetAIDAIon()->GetCluster(jj)->GetXEnergy();
                      aida.EY=evts->GetAIDAIon()->GetCluster(jj)->GetYEnergy();
                      tree->Fill();
                  }
              }
          }
          if ((Mode==4) && evts->IsBETA()){
              for (unsigned short jj=0;jj<evts->GetAIDABeta()->GetClusters().size();jj++){
                  if (evts->GetAIDABeta()->GetCluster(jj)->GetXMultiplicity()==1&&evts->GetAIDABeta()->GetCluster(jj)->GetYMultiplicity()==1){
                      aida.z=evts->GetAIDABeta()->GetCluster(jj)->GetHitPositionZ();
                      aida.y=evts->GetAIDABeta()->GetCluster(jj)->GetHitPositionY();
                      aida.x=evts->GetAIDABeta()->GetCluster(jj)->GetHitPositionX();
                      aida.EX=evts->GetAIDABeta()->GetCluster(jj)->GetXEnergy();
                      aida.EY=evts->GetAIDABeta()->GetCluster(jj)->GetYEnergy();
                      tree->Fill();
                  }
              }

          }


          //!Get run time
          if (evts->IsBETA()&&start==0) {
              if (evts->GetAIDABeta()->GetHits().size()>0) tstart = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
              start++;
          }else if (!evts->IsBETA()&&start==0){
              if (evts->GetAIDAIon()->GetHits().size()>0) tstart = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
              start++;
          }

          if(signal_received){
            break;
          }
      }

      if (evts->IsBETA()) {
          if (evts->GetAIDABeta()->GetHits().size()>0) tend = evts->GetAIDABeta()->GetHit(0)->GetTimestamp();
      }else if (!evts->IsBETA()){
          if (evts->GetAIDAIon()->GetHits().size()>0) tend = evts->GetAIDAIon()->GetHit(0)->GetTimestamp();
      }
      runtime[i+1] = (double)((tend-tstart)*ClockResolution)/(double)1e9;
      runtime[0] += runtime[i+1];
      cout<<evts->GetCurrentBetaEvent()<<" beta events"<<endl;
      cout<<evts->GetCurrentPulserEvent()<<" pulser events"<<endl;
      cout<<ttotal<<" all events"<<endl;

      delete evts;
  }

  if (Mode==3||Mode==4) {
      tree->Write();
  }
  else{
      spectra->Write();
      calibspectra->Write();
  }
  runtime.Write("runtime");

  ofile->Close();

  cout<<"\n**********************SUMMARY**********************\n"<<endl;
  cout<<"Total run length = "<<runtime[0]<< " seconds"<<endl;
  cout<<"Sub runs length"<<endl;
  for (Int_t i=0;i<nfiles;i++){
      cout<<inputfiles[i]<<" - "<<runtime[i+1]<< " seconds"<<endl;
  }
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
