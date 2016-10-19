#include "Riostream.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TMarker.h>
#include <TH2.h>
#include <TExec.h>
#include <TCanvas.h>

void combine_thr_y()
{
    Int_t dssd,strip;
    Double_t threshold[6][256];


    std::ifstream ifs("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            //clean theshold
            threshold[ds][i]=33000;
            ifs>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }
    cout<<"done step1"<<endl;
    std::ifstream ifsin("out_manual_threshold_final.txt");
    while (!ifsin.eof()){
        ifsin>>dssd>>strip>>threshold[dssd][strip];
        //cout<<threshold[dssd][strip]<<endl;
    }
    cout<<"done step2"<<endl;
    std::ifstream ifsin2("out_higher2_thr_y.txt");
    while (!ifsin2.eof()){
        ifsin2>>dssd>>strip>>threshold[dssd][strip];
        //cout<<threshold[dssd][strip]<<endl;
    }
    cout<<"done step3"<<endl;
    std::ofstream ofs("out_manual_threshold_final_higher2_y.txt");

    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            if (enable[ds][i]&&threshold[ds][i]==33000) cout<<"Missing strips:"<<ds<<" - "<<i<<endl;
            if (ds==5&&i==0) threshold[ds][i]=33000;//SPECIAL CASE
            ofs<<ds<<" "<<i<<" "<<threshold[ds][i]<<endl;
        }
    }

}


void manual_threshold_y()
{
    TFile *infile=new TFile("aida60_80_bkg_spec.root");
    TH2F* h2[6];
    for (Int_t i=0;i<6;i++)
    h2[i]=(TH2F*) infile->Get(Form("dssd%d",i));
    cout<<"Start over (y=yes,anykey=no)?"<<endl;
    TString ss;
    cin>>ss;
    if (ss=="y"){
        ofstream stemp("out.txt");
        stemp.close();
    }
    cout<<"Start dssd:"<<endl;
    Int_t startdssd;
    cin>>startdssd;
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;
    Int_t k=0;
    TH1F * histo[256];
    TCanvas* c1=new TCanvas("aa","aa",900,700);
    //Scan for available channels:
    std::ifstream ifs("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }
    Bool_t startFlag=false;
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=128;i<256;i++){
            if (ds==startdssd&&i==startCh) startFlag=true;
            if (enable[ds][i]&&startFlag){
                histo[i]=(TH1F*) h2[ds]->ProjectionY(Form("%d",i),i+1,i+1);
                histo[i]->GetXaxis()->SetRangeUser(0,1000);
                histo[i]->Draw();
                histo[i]->SetTitle(Form("%d",ds));
                c1->Update();
                TExec ex("ex1",".x exec1.C");
                ex.Draw();
                c1->WaitPrimitive();
                k++;
            }
        }
    }
    c1->WaitPrimitive();

}

void manual_threshold()
{
    TFile *infile=new TFile("aida60_80_bkg_spec.root");
    TH2F* h2[6];
    for (Int_t i=0;i<6;i++)
    h2[i]=(TH2F*) infile->Get(Form("dssd%d",i));
    cout<<"Start over (y=yes,anykey=no)?"<<endl;
    TString ss;
    cin>>ss;
    if (ss=="y"){
        ofstream stemp("out.txt");
        stemp.close();
    }
    cout<<"Start dssd:"<<endl;
    Int_t startdssd;
    cin>>startdssd;
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;
    Int_t k=0;
    TH1F * histo[256];
    TCanvas* c1=new TCanvas("aa","aa",900,700);
    //Scan for available channels:
    std::ifstream ifs("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }
    Bool_t startFlag=false;
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            if (ds==startdssd&&i==startCh) startFlag=true;
            if (enable[ds][i]&&startFlag){
                histo[i]=(TH1F*) h2[ds]->ProjectionY(Form("%d",i),i+1,i+1);
                histo[i]->Draw();
                histo[i]->SetTitle(Form("%d",ds));
                c1->Update();
                TExec ex("ex1",".x exec1.C");
                ex.Draw();
                c1->WaitPrimitive();
                k++;
            }
        }
    }
    c1->WaitPrimitive();

}


Double_t findThreshold(TH1F* spec,Int_t ncountsmax)
{
    Int_t ntotal=0;
    Double_t thr = 0;
    for (Int_t i=spec->GetNbinsX()+1;i>0;i--){
    //for (Int_t i=spec->GetNbinsX();i>0;i--){
        ntotal+=spec->GetBinContent(i);
        //if (i==spec->GetNbinsX()+1) cout<<ntotal<<endl;
        if (ntotal>ncountsmax) {
            thr = spec->GetBinCenter(i);
            break;
        }
    }
    return thr+0.5;
}

void auto_threshold(TString inputfile,TString outputfile, Int_t maxrate,Int_t sleeptime)
{

    TFile *infile=new TFile(inputfile);
    TH2F* h2[6];
    for (Int_t i=0;i<6;i++)
    h2[i]=(TH2F*) infile->Get(Form("dssd%d",i));
    TVectorD* runtime = (TVectorD*) infile->Get("runtime");
    Int_t ncountmax=(Int_t)maxrate* (Double_t)runtime(0);
    cout<<"Total run time: "<<(Double_t)runtime(0)<<" s"<<endl;
    cout<<"Start over (y=yes,anykey=no)?"<<endl;
    TString ss;
    cin>>ss;
    if (ss=="y"){
        ofstream stemp("out.txt");
        stemp.close();
    }
    cout<<"Start dssd:"<<endl;
    Int_t startdssd;
    cin>>startdssd;
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;
    Int_t k=0;
    TH1F * histo[256];
    TCanvas* c1=new TCanvas("aa","aa",900,700);

    //Scan for available channels:
    /*
    std::ifstream ifs("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }
    */
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
           enable[ds][i]=1;
        }
    }

    Bool_t startFlag=false;
    ofstream str(outputfile,ios::app);

    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            if (ds==startdssd&&i==startCh) startFlag=true;
            if (enable[ds][i]&&startFlag){
                histo[i]=(TH1F*) h2[ds]->ProjectionY(Form("%d",i),i+1,i+1);
                histo[i]->Draw();
                Double_t thr=findThreshold(histo[i],ncountmax);
		
                TMarker* mk=new TMarker(thr,histo[i]->GetBinContent(histo[i]->GetXaxis()->FindBin(thr)),1);
                mk->SetMarkerStyle(23);
                mk->SetMarkerSize(1.5);
                mk->SetMarkerColor(2);
                mk->Draw();
                if (thr>3000) thr=33000; //totally disable
		cout<<"thr = "<<thr<<endl;
                histo[i]->SetTitle(Form("%d",ds));
                c1->Update();
                str<<ds<<"\t"<<i<<"\t"<<"\t"<<thr<<endl;
                //str<<histo[i]->GetTitle()<<" "<<histo[i]->GetName()<<" "<<thr<<endl;
                gSystem->Sleep(sleeptime);
                k++;
            }
        }
    }
    str.close();
}


void refineThesholdTable()
{
    Int_t dssd,strip;
    Double_t threshold[6][256];


    std::ifstream ifs("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            //clean theshold
            threshold[ds][i]=33000;
            ifs>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }
    cout<<"done step1"<<endl;
    std::ifstream ifsin("out_manual_threshold.txt");
    while (!ifsin.eof()){
        ifsin>>dssd>>strip>>threshold[dssd][strip];
        //cout<<threshold[dssd][strip]<<endl;
    }

    std::ofstream ofs("out_manual_threshold_final.txt");

    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            if (enable[ds][i]&&threshold[ds][i]==33000) cout<<"Missing strips:"<<ds<<" - "<<i<<endl;
            ofs<<ds<<" "<<i<<" "<<threshold[ds][i]<<endl;
        }
    }
}
void refineFEETable()
{
    std::ifstream ifs2("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs2>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
        }
    }

    std::ifstream ifs("FEE_table_wdisabled.txt");
    Int_t DSSDtoFEE[6][256];
    Int_t DSSDtoCh[6][256];
    Int_t DSSDtoMask[6][256];
    Int_t fee,ch,dssd,strip,mask;
    while (!ifs.eof()){
        ifs>>fee>>ch>>dssd>>strip>>mask;
        DSSDtoFEE[dssd][strip]=fee;
        DSSDtoCh[dssd][strip]=ch;
        DSSDtoMask[dssd][strip]=mask;
    }

    std::ofstream ofs("FEE_table_refined.txt");
    Int_t k=0;
    for (Int_t i=0;i<6;i++){
        for (Int_t j=0;j<256;j++){
            ofs<<DSSDtoFEE[i][j]<<" "<<DSSDtoCh[i][j]<<"\t"<<i<<"\t"<<j<<"\t"<<enable[i][j]<<endl;
            if (enable[i][j]!=DSSDtoMask[i][j]&&DSSDtoMask[i][j]==0) {
                cout<<"Miss match: "<<i<<" - "<<j<<endl;
                k++;
            }

        }
    }
    cout<<"Total missmatched channel= "<<k<<endl;
    ofs.close();
}

void commonthresholdGen(TString outfile,Double_t keVx,Double_t keVy){
    Double_t threshold[6][256];
    std::ifstream ifs2("cal_table_final_safe.txt");
    Int_t dssd,ch;
    Double_t cal1,cal2;
    std::ofstream ofs(outfile);
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs2>>dssd>>ch>>cal1>>cal2;
            if (ch!=i||dssd!=ds) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) {
                threshold[ds][i]=33000;
            }else{
                if (i<128) threshold[ds][i]=(keVx-cal1)/cal2;
                else threshold[ds][i]=(keVy-cal1)/cal2;
                if (threshold[ds][i]<0) threshold[ds][i]=0;
            }
            ofs<<ds<<" "<<i<<" "<<threshold[ds][i]<<endl;
        }
    }
}

void getpeakwidth(TString outfile, Int_t msstop)
{
    std::ifstream ifs2("cal_table_final_safe.txt");
     std::ofstream str(outfile);
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable_d[6][256];
    Int_t enable[1536];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs2>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable_d[ds][i]=0;
            else enable_d[ds][i]=1;

            if (i<128) enable[ds*128+i]=enable_d[ds][i];
            else enable[768+ds*128+i-128]=enable_d[ds][i];
        }
    }

    TFile *f=TFile::Open("cal_specs63_safe_wcommonthr200_150keV.root");
    TH2F* sumspec=(TH2F*) f->Get("sumspec");
    TH1F* histo[1536];
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;

    Double_t sigma[1536];

    for (Int_t i=0;i<1536;i++) sigma[i]=33000/5; //5 sigma/no sense
    //special case
    sigma[250]=15.534;
    sigma[506]=15.94;
    sigma[588]=23.541;
    sigma[1000]=22.137;
    sigma[1001]=19.568;
    sigma[1219]=21.735;
    sigma[1267]=26.305;
    sigma[1269]=21.074;
    sigma[1270]=18.725;
    sigma[1273]=26.394;
    sigma[1274]=25.111;
    enable[250]=0;
    enable[506]=0;
    enable[588]=0;
    enable[1000]=0;
    enable[1001]=0;
    enable[1219]=0;
    enable[1267]=0;
    enable[1269]=0;
    enable[1270]=0;
    enable[1273]=0;
    enable[1274]=0;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    for (Int_t i=startCh;i<1536;i++){
        histo[i]=(TH1F*) sumspec->ProjectionY(Form("%d",i),i+1,i+1);
        if (enable[i]){
            histo[i]->Draw("colz");
            histo[i]->GetXaxis()->SetRangeUser(600,1200);
            histo[i]->Fit("gaus","QR","",600,1200);
            TF1* myfunc2=histo[i]->GetFunction("gaus");
            myfunc2->SetLineColor(3);
            myfunc2->Draw("LSAME");
            c1->Update();
            sigma[i]=myfunc2->GetParameter(2);
            cout<<sigma[i]<<endl;
            if (sigma[1536]*5>200) gSystem->Sleep(5000);
            gSystem->Sleep(msstop);
        }
    }
    for (Int_t i=0;i<6;i++){
        for (Int_t j=0;j<128;j++){
            Int_t chX=i*128+j;
            str<<i<<" "<<j<<" "<<sigma[chX]<<endl;
        }
        for (Int_t j=0;j<128;j++){
            Int_t chY=768+i*128+j;
            str<<i<<" "<<j+128<<" "<<sigma[chY]<<endl;
        }
    }
}

void make_sigma_threshold_table(TString sigmatable,TString caltable,TString outfile,Double_t factor)
{
    //read calibration table
    Double_t cal1[6][256];
    Double_t cal2[6][256];
    Double_t ccal1,ccal2;
    Int_t dssd,strip;
    std::ifstream ifs_cal(caltable);
    while (!ifs_cal.eof()){
        ifs_cal>>dssd>>strip>>ccal1>>ccal2;
        cal1[dssd][strip]=ccal1;
        cal2[dssd][strip]=ccal2;
    }

    std::ofstream str(outfile);
    std::ifstream ifs(sigmatable);
    Double_t sigma,threshold;
    while (!ifs.eof()){
        ifs>>dssd>>strip>>sigma;

        threshold=(sigma*factor-cal1[dssd][strip])/cal2[dssd][strip];
        if (threshold<0) threshold=0;
        if (sigma==33000/5) threshold=33000;
        str<<dssd<<" "<<strip<<" "<<threshold<<endl;
    }
}
