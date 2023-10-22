#include "ReaderTemplate.h"
#include "TH1.h"
void check1(TString infile){
    //! Init reader
    Reader* reader=new Reader;
    reader->InitBeta(infile);

    //! book histograms and trees
    TH1F* histo=new TH1F("mult","mult",10,0,10);

    //! get number of entries
    Long64_t nentries = reader->GetEntries();
    AIDA* aida = reader->aida;
    //! event loop
    for (Long64_t jentry = 0;jentry<nentries;jentry++){
        reader->GetEntry(jentry);
        histo->Fill(aida->GetMult());
    }
    histo->Draw("colz");
}
