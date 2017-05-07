
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
   double E;
   uint64_t T;
   uint16_t Id;
   uint16_t type;
   uint16_t Index1;
   uint16_t Index2;
   uint16_t InfoFlag;
   std::string Name;
   void clear();
};

#endif
