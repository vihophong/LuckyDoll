#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <TString.h>
#include <string>
#include <TSpectrum.h>
#include <TText.h>
#include <TF1.h>
#include "TExec.h"
#include <vector>
#include "TSystem.h"
#include "TGraph.h"
#define NumDSSD 6
#define NumStrX 128
#define NumStrY 128
#define NumStrXY 256
#define NumFee 33
#define NumChFee 64

void makespecs(TString outfile){
    ifstream ifs("file_list.txt");
    Int_t nline;
    ifs>>nline;
    TFile *f;
    TString filelist[nline];
    TH2F* specs[100];
    TH2F* sumSpec;
    sumSpec=new TH2F ("spectra","spectra",1536,0,1536,4000,0,24000);
    for (Int_t i=0;i<nline;i++){
        ifs>>filelist[i];
        //f=TFile::Open("aidasecond/"+filelist[i]);
        f=TFile::Open("bicalsecond/"+filelist[i]);
        specs[i]=(TH2F*) f->Get("spectra");
        specs[i]->SetName(Form("spec%d",i));
        specs[i]->SetTitle(filelist[i]);
        sumSpec->Add(specs[i]);
    }
    sumSpec->SetName("sumspec");
    sumSpec->Draw("colz");
    TFile* outspecs=new TFile(outfile,"recreate");
    sumSpec->Write();
    for (Int_t i=0;i<nline;i++) specs[i]->Write();
    outspecs->Close();
}
void viewbadspec(TString badlist){
    std::ifstream ifs(badlist);
    Int_t chnum;
    TCanvas* spec;
    while (!ifs.eof()){
        ifs>>chnum;
        viewSpec("specs63.root",chnum);
        gSystem->Sleep(5000);
    }
}
void viewSpeclowE(TString inputfile,Int_t sch,Double_t standard=0.01)
{
    Int_t libpeaks[9]={90000,80000,70000,60000,50000,40000,30000,20000,10000};
    TFile *f=TFile::Open(inputfile);
    TH2F* sumspec=(TH2F*) f->Get("sumspec");
    TH1F* histo[1536];
    Int_t k=0;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    Int_t i=sch;
    histo[i]=(TH1F*)sumspec->ProjectionY(Form("%d",i),i+1,i+1);
    //list all peaks
    TSpectrum *s=new TSpectrum();
    s->Search(histo[i],10,"",standard);
    Float_t* x=s->GetPositionX();
    Int_t NPeaks=s->GetNPeaks();
    std::vector<float> xvec;
    for (Int_t j=0;j<NPeaks;j++) xvec.push_back(x[j]);
    std::sort (xvec.begin(),xvec.end());
    //get only first 9 peaks
    Double_t valpeak[9];
    Int_t n=0;
    for (Int_t j=NPeaks-1;j>=0;j--){
        if (j>=0&&n<4) {
            //if ((xvec[j]-xvec[j-1])>500){
                Double_t val = histo[i]->GetBinContent(histo[i]->GetXaxis()->FindBin(xvec[j]));
                TText *t = new TText(xvec[j],val,Form("%d",libpeaks[n]));
                t->SetTextSize(0.02);
                t->Draw();
                cout<<n<<"-"<<xvec[j]<<"-"<<val<<endl;
                valpeak[n]=xvec[j];
                n++;
            //}
        }
        if (n>=4) break;
    }
    for (Int_t i=0;i<4;i++) cout<<valpeak[i]<<" ";
    cout<<endl;
    c1->Update();
    /*
    TExec ex("ex1",".x exec1.C");
    ex.Draw();
    c1->WaitPrimitive();
    */
    k++;
    //TSystem::Sleep(1000);
    //gSystem->Sleep(1000);

}

void viewSpec(TString inputfile,Int_t sch,Double_t standard=0.01)
{
    Int_t libpeaks[4]={10000,8000,6000,4000};
    TFile *f=TFile::Open(inputfile);
    TH2F* sumspec=(TH2F*) f->Get("sumspec");
    TH1F* histo[1536];
    Int_t k=0;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    Int_t i=sch;
    histo[i]=(TH1F*)sumspec->ProjectionY(Form("%d",i),i+1,i+1);
    //list all peaks
    TSpectrum *s=new TSpectrum();
    s->Search(histo[i],20,"",standard);
    Float_t* x=s->GetPositionX();
    Int_t NPeaks=s->GetNPeaks();
    std::vector<float> xvec;
    for (Int_t j=0;j<NPeaks;j++) xvec.push_back(x[j]);
    std::sort (xvec.begin(),xvec.end());
    //get only first 9 peaks
    Double_t valpeak[9];
    Int_t n=0;
    for (Int_t j=NPeaks-1;j>=0;j--){
        if (j>0) {
            //if ((xvec[j]-xvec[j-1])>500){
                histo[i]->GetXaxis()->SetRangeUser(0,5000);
                Double_t val = histo[i]->GetBinContent(histo[i]->GetXaxis()->FindBin(xvec[j]));
                TText *t = new TText(xvec[j],val,Form("%d",libpeaks[n]));
                t->SetTextSize(0.02);
                t->Draw();
                cout<<n<<"-"<<xvec[j]<<"-"<<val<<endl;
                valpeak[n]=xvec[j];
                n++;
            //}
        }else if (j==0&&n<9){
            cout<<"sss"<<endl;
        }
        if (n>=4) break;
    }
    c1->Update();
    /*
    TExec ex("ex1",".x exec1.C");
    ex.Draw();
    c1->WaitPrimitive();
    */
    k++;
    //TSystem::Sleep(1000);
    //gSystem->Sleep(1000);

}
void makePeaksTable_lowE(TString inputfile)
{
    std::ifstream ifs2("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;
    Int_t enable[6][256];
    Int_t enablel[1536];
    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs2>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            if (cal1==0&&cal2==0) enable[ds][i]=0;
            else enable[ds][i]=1;
            Int_t index;
            if (i<128) index=ds*128+i;
            else index=768+ds*128+i-128;
            enablel[index]=enable[ds][i];
        }
    }

    cout<<"Start over (y=yes,anykey=no)?"<<endl;
    TString ss;
    cin>>ss;
    if (ss=="y"){
        ofstream stemp("out_lowE.txt");
        stemp.close();
        ofstream stemp2("out_lowE_3.txt");
        stemp2.close();
    }
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;

    ofstream str("out_lowE.txt",ios::app);
    ofstream str_3("out_lowE_3.txt",ios::app);
    std::ofstream ofs_bad("out_bad_peaks_lowE.txt");

    Int_t libpeaks[4]={10000,8000,6000,4000};
    TFile *f=TFile::Open(inputfile);
    TH2F* sumspec=(TH2F*) f->Get("sumspec");
    TFile *outspecs=new TFile("out.root","recreate");
    TH1F* histo[1536];
    Int_t k=0;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    for (Int_t i=startCh;i<1536;i++){
        if (!enablel[i]) continue;

        histo[i]=(TH1F*) sumspec->ProjectionY(Form("%d",i),i+1,i+1);
        //list all peaks
        histo[i]->Draw("hist");
        histo[i]->GetXaxis()->SetRangeUser(0,5000);
        TSpectrum *s=new TSpectrum();
        s->Search(histo[i],10,"",0.005);
        Float_t* x=s->GetPositionX();
        Int_t NPeaks=s->GetNPeaks();
        std::vector<float> xvec;
        for (Int_t j=0;j<NPeaks;j++) xvec.push_back(x[j]);
        std::sort (xvec.begin(),xvec.end());
        //get only first 4 peaks
        Double_t valpeak[4];
        Int_t n=0;
        for (Int_t j=NPeaks-1;j>=0;j--){
            if (j>=0&&n<4) {
                //if ((xvec[j]-xvec[j-1])>50){
                    Double_t val = histo[i]->GetBinContent(histo[i]->GetXaxis()->FindBin(xvec[j]));
                    TText *t = new TText(xvec[j],val,Form("%d",libpeaks[n]));
                    t->SetTextSize(0.02);
                    t->Draw();
                    cout<<NPeaks<<"-"<<n<<"-"<<xvec[j]<<"-"<<val<<endl;
                    valpeak[n]=xvec[j];
                    n++;
                //}
            }
            //else if (j==0&&n<3){
            //    cout<<"sss"<<i<<endl;
            //}
            if (n>=4) break;
        }
        if (n==4){//found 4 peaks
            str<<i<<" ";
            for (Int_t j=0;j<4;j++){
                str<<valpeak[j]<<" ";
            }
            str<<endl;
        }else if (n==3){//only found 3 peaks
            str_3<<i<<" ";
            for (Int_t j=0;j<3;j++){
                str_3<<valpeak[j]<<" ";
            }
            str_3<<endl;
        }else{
            ofs_bad<<i<<endl;
        }
        c1->Update();
        c1->SetName(Form("his%d",i));
        c1->Write();
        /*
        TExec ex("ex1",".x exec1.C");
        ex.Draw();
        c1->WaitPrimitive();
        */
        k++;
        //TSystem::Sleep(1000);
        gSystem->Sleep(3000);
    }
    //c1->WaitPrimitive();
    outspecs->Close();
}
void makePeaksTable(TString inputfile)
{
    cout<<"Start over (y=yes,anykey=no)?"<<endl;
    TString ss;
    cin>>ss;
    if (ss=="y"){
        ofstream stemp("out.txt");
        stemp.close();
    }
    cout<<"Start channel:"<<endl;
    Int_t startCh;
    cin>>startCh;

    ofstream str("out.txt",ios::app);
    std::ofstream ofs_bad("out_bad_peaks.txt");

    Int_t libpeaks[9]={90000,80000,70000,60000,50000,40000,30000,20000,10000};
    TFile *f=TFile::Open(inputfile);
    TH2F* sumspec=(TH2F*) f->Get("sumspec");
    TFile *outspecs=new TFile("out.root","recreate");
    TH1F* histo[1536];
    Int_t k=0;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    for (Int_t i=startCh;i<1536;i++){
        histo[i]=(TH1F*) sumspec->ProjectionY(Form("%d",i),i+1,i+1);
        //list all peaks
        TSpectrum *s=new TSpectrum();
        s->Search(histo[i],20,"",0.01);
        Float_t* x=s->GetPositionX();
        Int_t NPeaks=s->GetNPeaks();
        std::vector<float> xvec;
        for (Int_t j=0;j<NPeaks;j++) xvec.push_back(x[j]);
        std::sort (xvec.begin(),xvec.end());
        //get only first 9 peaks
        Double_t valpeak[9];
        Int_t n=0;
        for (Int_t j=NPeaks-1;j>=0;j--){
            if (j>0) {
                if ((xvec[j]-xvec[j-1])>500){
                    Double_t val = histo[i]->GetBinContent(histo[i]->GetXaxis()->FindBin(xvec[j]));
                    TText *t = new TText(xvec[j],val,Form("%d",libpeaks[n]));
                    t->SetTextSize(0.02);
                    t->Draw();
                    cout<<n<<"-"<<xvec[j]<<"-"<<val<<endl;
                    valpeak[n]=xvec[j];
                    n++;
                }
            }else if (j==0&&n<9){
                cout<<"sss"<<endl;
            }
            if (n>=9) break;
        }
        if (n==9){
            str<<i<<" ";
            for (Int_t j=0;j<9;j++){
                str<<valpeak[j]<<" ";
            }
            str<<endl;
        }else{
            ofs_bad<<i<<endl;
        }
        c1->Update();
        c1->SetName(Form("his%d",i));
        c1->Write();
        /*
        TExec ex("ex1",".x exec1.C");
        ex.Draw();
        c1->WaitPrimitive();
        */
        k++;
        //TSystem::Sleep(1000);
        //gSystem->Sleep(1000);
    }
    //c1->WaitPrimitive();
    outspecs->Close();
}
void calcurvelowE(TString infile3,TString infile4,TString outfile){
    std::ofstream str(outfile);
    std::ifstream ifs3(infile3);
    std::ifstream ifs(infile4);
    std::ofstream ofs_bad("bad_cal_list.txt");
    Int_t k=0;
    Int_t kk=0;
    Double_t libpeaks[4]={10000,8000,6000,4000};
    Double_t peakval[4];
    Double_t Ecal_pol1[5];
    Double_t Ecal2_pol1[5];

    Int_t ch;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    TH1F* ee=new TH1F("ss","ss",5000,0.99,1.);
    Double_t vperch[1536];
    Double_t choffset[1536];
    Double_t kevperch[1536];
    Double_t offset[1536];
    Double_t rsquare[1536];



    for (Int_t i=0;i<1536;i++){
        vperch[i]=0;
        choffset[i]=0;
        kevperch[i]=0;
        offset[i]=0;
    }


    std::ifstream ifs2("cal_table_final_safe.txt");
    Int_t temp;
    Double_t cal1,cal2;

    for (Int_t ds=0;ds<6;ds++){
        for (Int_t i=0;i<256;i++){
            ifs2>>temp>>temp>>cal1>>cal2;
            if (temp!=i) cout<<"Warning at "<<ds<<"-"<<i<<endl;
            Int_t index;
            if (i<128) index=ds*128+i;
            else index=768+ds*128+i-128;
            kevperch[index]=cal2;
            offset[index]=cal1;
        }
    }



    //3 peaks
    while (!ifs3.eof()){
        ifs3>>ch;
        for (Int_t i=0;i<3;i++) ifs3>>peakval[i];
        TGraph *g=new TGraph(3,peakval,libpeaks);
        TGraph *grv=new TGraph(3,libpeaks,peakval);
        g->Draw("AP*");
        g->Fit("pol1","Q");
        grv->Fit("pol1","Q");
        TF1* myfunc2=g->GetFunction("pol1");
        myfunc2->SetLineColor(3);
        myfunc2->Draw("LSAME");
        Ecal_pol1[0]=myfunc2->GetParameter(0);
        Ecal_pol1[1]=myfunc2->GetParameter(1);
        Ecal_pol1[2]=myfunc2->GetParError(0);
        Ecal_pol1[3]=myfunc2->GetParError(1);
        Ecal_pol1[4]=myfunc2->GetChisquare()/myfunc2->GetNDF();

        TF1* myfunc3=grv->GetFunction("pol1");
        Ecal2_pol1[0]=myfunc3->GetParameter(0);
        Ecal2_pol1[1]=myfunc3->GetParameter(1);
        Ecal2_pol1[2]=myfunc3->GetParError(0);
        Ecal2_pol1[3]=myfunc3->GetParError(1);
        Ecal2_pol1[4]=myfunc3->GetChisquare()/myfunc3->GetNDF();

        //Calculate regression sum of square
        Double_t sum1=0;
        Double_t sum2=0;
        for (Int_t i=0;i<3;i++){
            sum1+=(libpeaks[i]-50000)*(libpeaks[i]-50000);
            sum2+=(libpeaks[i]-myfunc2->Eval(peakval[i]))*(libpeaks[i]-myfunc2->Eval(peakval[i]));
        }
        rsquare[ch]=(sum1-sum2)/sum1;

        TText *t = new TText(peakval[2],libpeaks[2],Form("%f",rsquare[ch]));
        t->SetTextSize(0.05);
        t->Draw();
        TText *tt = new TText(peakval[1],libpeaks[1],Form("ch= %d",ch));
        tt->SetTextSize(0.05);
        tt->Draw();
        c1->Update();
        k++;
        //if (rsquare[ch]<0.999) {//safe
        //if (rsquare[ch]<0.96) { //dangerous
        if (rsquare[ch]<0.96) { //next level of dangerous
            //if (Ecal_pol1[4]>300000) {
            ofs_bad<<ch<<" "<<rsquare[ch]<<endl;
            gSystem->Sleep(10000);
            kk++;
        }else{
            choffset[ch]=Ecal2_pol1[0];
            vperch[ch]=Ecal_pol1[1];
        }
        ee->Fill(rsquare[ch]);
    }
    cout<<"finished 3 peaks calib"<<endl;
    //4 peaks
    while (!ifs.eof()){
        ifs>>ch;
        for (Int_t i=0;i<4;i++) ifs>>peakval[i];
        TGraph *g=new TGraph(4,peakval,libpeaks);
        TGraph *grv=new TGraph(4,libpeaks,peakval);
        g->Draw("AP*");
        g->Fit("pol1","Q");
        grv->Fit("pol1","Q");
        TF1* myfunc2=g->GetFunction("pol1");
        myfunc2->SetLineColor(3);
        myfunc2->Draw("LSAME");
        Ecal_pol1[0]=myfunc2->GetParameter(0);
        Ecal_pol1[1]=myfunc2->GetParameter(1);
        Ecal_pol1[2]=myfunc2->GetParError(0);
        Ecal_pol1[3]=myfunc2->GetParError(1);
        Ecal_pol1[4]=myfunc2->GetChisquare()/myfunc2->GetNDF();

        TF1* myfunc3=grv->GetFunction("pol1");
        Ecal2_pol1[0]=myfunc3->GetParameter(0);
        Ecal2_pol1[1]=myfunc3->GetParameter(1);
        Ecal2_pol1[2]=myfunc3->GetParError(0);
        Ecal2_pol1[3]=myfunc3->GetParError(1);
        Ecal2_pol1[4]=myfunc3->GetChisquare()/myfunc3->GetNDF();

        //Calculate regression sum of square
        Double_t sum1=0;
        Double_t sum2=0;
        for (Int_t i=0;i<4;i++){
            sum1+=(libpeaks[i]-50000)*(libpeaks[i]-50000);
            sum2+=(libpeaks[i]-myfunc2->Eval(peakval[i]))*(libpeaks[i]-myfunc2->Eval(peakval[i]));
        }
        rsquare[ch]=(sum1-sum2)/sum1;

        TText *t = new TText(peakval[3],libpeaks[3],Form("%f",rsquare[ch]));
        t->SetTextSize(0.05);
        t->Draw();
        TText *tt = new TText(peakval[2],libpeaks[2],Form("ch= %d",ch));
        tt->SetTextSize(0.05);
        tt->Draw();
        c1->Update();
        k++;
        //if (rsquare[ch]<0.999) {//safe
        //if (rsquare[ch]<0.96) { //dangerous
        if (rsquare[ch]<0.96) { //next level of dangerous
            //if (Ecal_pol1[4]>300000) {
            ofs_bad<<ch<<" "<<rsquare[ch]<<endl;
            gSystem->Sleep(10000);
            kk++;
        }else{
            choffset[ch]=Ecal2_pol1[0];
            vperch[ch]=Ecal_pol1[1];
        }
        ee->Fill(rsquare[ch]);
    }

    //get average V/ch
    Double_t totalvperchX,totalvperchY;
    totalvperchX=totalvperchY=0;
    Int_t nchX,nchY;
    nchX=nchY=0;
    for (Int_t i=0;i<1536;i++){
        if (!(vperch[i]==0&&choffset[i]==0)){
            if (i<768){
                totalvperchX+=vperch[i];
                nchX++;
            }else{
                totalvperchY+=vperch[i];
                nchY++;
            }
        }
    }
    Double_t averagevperchX=totalvperchX/nchX;
    Double_t averagevperchY=totalvperchY/nchY;
    Double_t xecal=0.753271;
    cout<<"A:"<<totalvperchX<<"-"<<totalvperchY<<"-"<<nchX<<"-"<<nchY<<"-/"<<averagevperchX<<"-"<<averagevperchY<<endl;
    //no manual cal

    //scan for manual calibration (before this step)
    std::ifstream ifs_manual("manual_callowE.txt");
    while (!ifs_manual.eof()){
        ifs_manual>>ch;
        ifs_manual>>vperch[ch]>>choffset[ch];
        choffset[ch]=-choffset[ch]/vperch[ch];
        cout<<"manual: "<<ch<<"-"<<vperch[ch]<<"-"<<choffset[ch]<<endl;
    }


    Int_t ngood=0;
    Int_t nbad=0;
    std::ofstream ofs_nospec("zeroCh.txt");
    //reject channel with poor resolution
    Int_t rejectlist[11]={120,328,385,661,831,845,961,1012,1140,1229,1396};
    for (Int_t i=0;i<11;i++){
        vperch[rejectlist[i]]=0;choffset[rejectlist[i]]=0;
    }

    for (Int_t i=0;i<1536;i++){
        if (!(vperch[i]==0&&choffset[i]==0)){
            if (i<768) kevperch[i]=vperch[i]*xecal/averagevperchX;
            else kevperch[i]=vperch[i]*xecal/averagevperchX;
            offset[i]=-choffset[i]*kevperch[i];
            ngood++;
            //cout<<i<<"-"<<kevperch[i]<<"-"<<offset[i]<<endl;
        }else{
            ofs_nospec<<i<<endl;
            nbad++;
        }
    }
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrX;j++){
            Int_t chX=i*NumStrX+j;
            str<<i<<" "<<j<<" "<<offset[chX]<<" "<<kevperch[chX]<<endl;
        }
        for (Int_t j=0;j<NumStrY;j++){
            Int_t chY=NumDSSD*NumStrX+i*NumStrY+j;
            str<<i<<" "<<j+NumStrX<<" "<<offset[chY]<<" "<<kevperch[chY]<<endl;
        }
    }
    //str<<i/256<<" "<<i%256<<" "<<offset[i]<<" "<<kevperch[i]<<endl;

    cout<<"processs "<<ngood<<" spec"<<endl;
    cout<<"there are "<<nbad<<" bad spec"<<endl;
    ee->Draw("hist");

}
void calcurve(TString infile,TString outfile){
    std::ofstream str(outfile);
    std::ifstream ifs(infile);
    std::ofstream ofs_bad("bad_cal_list.txt");
    Int_t k=0;
    Int_t kk=0;
    Double_t libpeaks[9]={90000,80000,70000,60000,50000,40000,30000,20000,10000};
    Double_t peakval[9];
    Double_t Ecal_pol1[5];
    Double_t Ecal2_pol1[5];

    Int_t ch;
    TCanvas* c1=new TCanvas("hisname","hisname",800,600);
    TH1F* ee=new TH1F("ss","ss",5000,0.99,1.);
    Double_t vperch[1536];
    Double_t choffset[1536];
    Double_t kevperch[1536];
    Double_t offset[1536];
    Double_t rsquare[1536];
    for (Int_t i=0;i<1536;i++){
        vperch[i]=0;
        choffset[i]=0;
        kevperch[i]=0;
        offset[i]=0;
    }
    while (!ifs.eof()){
        ifs>>ch;
        for (Int_t i=0;i<9;i++) ifs>>peakval[i];
        TGraph *g=new TGraph(9,peakval,libpeaks);
        TGraph *grv=new TGraph(9,libpeaks,peakval);
        g->Draw("AP*");
        g->Fit("pol1","Q");
        grv->Fit("pol1","Q");
        TF1* myfunc2=g->GetFunction("pol1");
        myfunc2->SetLineColor(3);
        myfunc2->Draw("LSAME");
        Ecal_pol1[0]=myfunc2->GetParameter(0);
        Ecal_pol1[1]=myfunc2->GetParameter(1);
        Ecal_pol1[2]=myfunc2->GetParError(0);
        Ecal_pol1[3]=myfunc2->GetParError(1);
        Ecal_pol1[4]=myfunc2->GetChisquare()/myfunc2->GetNDF();

        TF1* myfunc3=grv->GetFunction("pol1");
        Ecal2_pol1[0]=myfunc3->GetParameter(0);
        Ecal2_pol1[1]=myfunc3->GetParameter(1);
        Ecal2_pol1[2]=myfunc3->GetParError(0);
        Ecal2_pol1[3]=myfunc3->GetParError(1);
        Ecal2_pol1[4]=myfunc3->GetChisquare()/myfunc3->GetNDF();

        //Calculate regression sum of square
        Double_t sum1=0;
        Double_t sum2=0;
        for (Int_t i=0;i<9;i++){
            sum1+=(libpeaks[i]-50000)*(libpeaks[i]-50000);
            sum2+=(libpeaks[i]-myfunc2->Eval(peakval[i]))*(libpeaks[i]-myfunc2->Eval(peakval[i]));
        }
        rsquare[ch]=(sum1-sum2)/sum1;

        TText *t = new TText(peakval[8],libpeaks[8],Form("%f",rsquare[ch]));
        t->SetTextSize(0.05);
        t->Draw();
        TText *tt = new TText(peakval[7],libpeaks[7],Form("ch= %d",ch));
        tt->SetTextSize(0.05);
        tt->Draw();
        c1->Update();
        k++;
        //if (rsquare[ch]<0.999) {//safe
        if (rsquare[ch]<0.96) { //notsafe
            //if (Ecal_pol1[4]>300000) {
            ofs_bad<<ch<<" "<<rsquare[ch]<<endl;
            //gSystem->Sleep(5000);
            kk++;
        }else{
            choffset[ch]=Ecal2_pol1[0];
            vperch[ch]=Ecal_pol1[1];
        }
        ee->Fill(rsquare[ch]);
    }
    //get average V/ch
    Double_t totalvperchX,totalvperchY;
    totalvperchX=totalvperchY=0;
    Int_t nchX,nchY;
    nchX=nchY=0;
    for (Int_t i=0;i<1536;i++){
        if (!(vperch[i]==0&&choffset[i]==0)){
            if (i<768){
                totalvperchX+=vperch[i];
                nchX++;
            }else{
                totalvperchY+=vperch[i];
                nchY++;
            }
        }
    }
    Double_t averagevperchX=totalvperchX/nchX;
    Double_t averagevperchY=totalvperchY/nchY;
    Double_t xecal=0.753271;
    cout<<"A:"<<totalvperchX<<"-"<<totalvperchY<<"-"<<nchX<<"-"<<nchY<<"-/"<<averagevperchX<<"-"<<averagevperchY<<endl;
    //scan for manual calibration (before this step)
    std::ifstream ifs_manual("manual_cal.txt");
    while (!ifs_manual.eof()){
        ifs_manual>>ch;
        ifs_manual>>vperch[ch]>>choffset[ch];
        choffset[ch]=-choffset[ch]/vperch[ch];
        cout<<"manual: "<<ch<<"-"<<vperch[ch]<<"-"<<choffset[ch]<<endl;
    }

    Int_t ngood=0;
    Int_t nbad=0;
    std::ofstream ofs_nospec("zeroCh.txt");
    //reject channel with poor resolution
    Int_t rejectlist[11]={120,328,385,661,831,845,961,1012,1140,1229,1396};
    for (Int_t i=0;i<11;i++){
        vperch[rejectlist[i]]=0;choffset[rejectlist[i]]=0;
    }

    for (Int_t i=0;i<1536;i++){
        if (!(vperch[i]==0&&choffset[i]==0)){
            if (i<768) kevperch[i]=vperch[i]*xecal/averagevperchX;
            else kevperch[i]=vperch[i]*xecal/averagevperchX;
            offset[i]=-choffset[i]*kevperch[i];
            ngood++;
            //cout<<i<<"-"<<kevperch[i]<<"-"<<offset[i]<<endl;
        }else{
            ofs_nospec<<i<<endl;
            nbad++;
        }
    }
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrX;j++){
            Int_t chX=i*NumStrX+j;
            str<<i<<" "<<j<<" "<<offset[chX]<<" "<<kevperch[chX]<<endl;
        }
        for (Int_t j=0;j<NumStrY;j++){
            Int_t chY=NumDSSD*NumStrX+i*NumStrY+j;
            str<<i<<" "<<j+NumStrX<<" "<<offset[chY]<<" "<<kevperch[chY]<<endl;
        }
    }
    //str<<i/256<<" "<<i%256<<" "<<offset[i]<<" "<<kevperch[i]<<endl;

    cout<<"processs "<<ngood<<" spec"<<endl;
    cout<<"there are "<<nbad<<" bad spec"<<endl;
    ee->Draw("hist");

}

void makeCalTable(TString infile,TString outfile){
    //std::ofstream str(outfile);
    std::ifstream ifs(infile);
    //first loop determine average kev/ch
    Int_t ch,nch;
    nch=0;
    Double_t temp;
    Double_t kevperch;
    Double_t totalkevperch;
    while (!ifs.eof()){
        ifs>>ch>>temp>>kevperch>>temp>>temp>>temp>>temp>>temp>>temp>>temp>>temp;
        cout<<ch<<"-"<<kevperch<<endl;
        totalkevperch+=kevperch;
        nch++;
    }
    Double_t average=totalkevperch/nch;
    cout<<"average="<<average<<endl;
    ifs.clear();
    ifs.seekg(0);
    Double_t xecal=0.753271;
    Double_t yecal=0.795423;
    Double_t offset;


}
