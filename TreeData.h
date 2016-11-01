
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>

#ifndef TreeData_h
#define TreeData_h
class TreeData  {

public:
   TreeData() { }
   ~TreeData() {}
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
