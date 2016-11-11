
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>

#ifndef BrikenTreeData_h
#define BrikenTreeData_h
class BrikenTreeData  {

public:
   BrikenTreeData() { }
   ~BrikenTreeData() {}
   Double_t        E;
   ULong64_t       T;
   UShort_t        Id;
   UShort_t        type;
   UShort_t        Index1;
   UShort_t        Index2;
   UShort_t        InfoFlag;
   std::string          Name;
};

#endif
