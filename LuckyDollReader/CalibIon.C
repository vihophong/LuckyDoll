#include "ReaderTemplate.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
void CalibIon(TString infile){
    Reader* reader=new Reader;
    reader->InitBeta(infile);

    TH2F* hisxy[NumDSSD];
    for (Int_t i=0;i<NumDSSD;i++) hisxy[i]=new TH2F(Form("hisXY-XY-%d",i),Form("hisXY-XY-%d",i),256,0,256,256,0,256);

    Long64_t nentries = reader->GetEntries();
    for (Long64_t jentry = 0;jentry<nentries;jentry++){
        reader->GetEntry(jentry);
        Double_t percent_complete=(Double_t)jentry/(Double_t)nentries*100;
        AIDA* ad = reader->aida;
        if (jentry%3000==0) cout<<percent_complete<<" % data piece proceeded"<<endl;

        Int_t posX[NumDSSD][5];
        Double_t EX[NumDSSD][5];
        Double_t maxEX[NumDSSD];
        Int_t maxEposX[NumDSSD];
        Int_t multX[NumDSSD];

        Int_t posY[NumDSSD][5];
        Double_t EY[NumDSSD][5];
        Double_t maxEY[NumDSSD];
        Int_t maxEposY[NumDSSD];
        Int_t multY[NumDSSD];

        for (int i=0;i<NumDSSD;i++){
            maxEX[i]=-9999;
            maxEposX[i]=-9999;
            multX[i]=0;
            maxEY[i]=-9999;
            maxEposY[i]=-9999;
            multY[i]=0;
        }


        for (int i=0;i<ad->GetMult();i++){
            AIDAHit* hit = ad->GetHit(i);
            Int_t z=hit->GetZ();
            //if (hit->GetZ()==3){

                //! x strips
                if (hit->GetXY()<128){//select x
                    Double_t ex=(Double_t) hit->GetEnergy();
                    if (ex>0){
                        if (multX[z]<5){
                            posX[z][multX[z]] = hit->GetXY();
                            EX[z][multX[z]] = ex;
                            if (ex>maxEX[z]){
                                maxEX[z]=ex;
                                maxEposX[z]=hit->GetXY();
                            }
                        }
                        multX[z]++;
                    }
                }
                //! y strips
                if (hit->GetXY()>128){//select x
                    Double_t ey=(Double_t) hit->GetEnergy();
                    if (ey>0){
                        if (multY[z]<5){
                            posY[z][multY[z]] = hit->GetXY();
                            EY[z][multY[z]] = ey;
                            if (ey>maxEY[z]){
                                maxEY[z]=ey;
                                maxEposY[z]=hit->GetXY();
                            }
                        }
                        multY[z]++;
                    }
                }
            //}//select dssd 3
        }
        for (Int_t z=0;z<NumDSSD;z++){
            if (multX[z]<5&&multX[z]>1){ //only select multipicity that larger than 2 in X
                for (Int_t j=0;j<multX[z];j++){
                    if (posX[z][j]!=maxEposX[z]){
                        hisxy[z]->Fill(maxEposX[z],posX[z][j]);
                    }
                }
            }
            if (multY[z]<5&&multY[z]>1){ //only select multipicity that larger than 2 in X
                for (Int_t j=0;j<multY[z];j++){
                    if (posY[z][j]!=maxEposY[z]){
                        hisxy[z]->Fill(maxEposY[z],posY[z][j]);
                    }
                }
            }
        }
        if (percent_complete>10) break;
    }
    TCanvas* c1=new TCanvas("crosstalk","crosstalk",900,700);
    c1->Divide(3,2);
    for (Int_t i=0;i<NumDSSD;i++){
        c1->cd(i+1);
        hisxy[i]->Draw("colz");
    }
}

