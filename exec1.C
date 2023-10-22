#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "Riostream.h"
#include <TROOT.h>
#include "TExec.h"
void exec1()
{
  //example of macro called when a pad is redrawn
  //one must create a TExec object in the following way
  // TExec ex("ex",".x exec1.C");
  // ex.Draw();
  // this macro prints the bin number and the bin content when one clicks
  //on the histogram contour of any histogram in a pad
  //Author: Rene Brun
  ofstream str("reject_list.txt",ios::app);
  int event = gPad->GetEvent();
  if (event != 11) return;
  int px = gPad->GetEventX();
  TObject *select = gPad->GetSelected();
  if (!select) return;
  if (select->InheritsFrom("TH1")) {
    TH1 *h = (TH1*)select;
    Float_t xx = gPad->AbsPixeltoX(px);
    Float_t x  = gPad->PadtoX(xx);
    Int_t binx = h->GetXaxis()->FindBin(x);
    str<<h->GetName<<endl;
    printf("event=%d, hist:%s, x=%f, content=%f\n",event,h->GetName(),x,h->GetBinContent(binx));
  }
  str.close();
}
