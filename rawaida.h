//////////////////////////////////////////////////////////
// Modified from old version of aida2tree
// Coded by Vi Ho Phong - phong@ribf.riken.jp
// Last modified: Sep 26,2016
//////////////////////////////////////////////////////////

#ifndef rawaida_h
#define rawaida_h

#include "Riostream.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>

class rawmidas_info{
public:
    rawmidas_info(){}
    virtual ~rawmidas_info(){}

    //in principal the data should not exceed 32 bit = 4byte =
    unsigned int timestampLsb;  //least significant bits timestamp
    unsigned int infoField;

    //  unsigned short sample_data[4];
    unsigned int adcData;
    unsigned int dataType;
    // type: 0= sample waveform, 1= sample lenght, 2= info data, 3= ADC data
    unsigned int feeId;
    unsigned int chId;
    unsigned int adcRange;
    unsigned int infoCode;
    //sample lenght ->omit
    unsigned int sampleLengh;
    unsigned int waveformtrace1;
    unsigned int waveformtrace2;
    unsigned int waveformtrace3;
    unsigned int waveformtrace4;
    void  Show(){
        cout<<"timestampLsb "<<timestampLsb<<endl;
        cout<<"infoField "<<infoField<<endl;
        cout<<"adcData "<<adcData<<endl;
        cout<<"dataType "<<dataType<<endl;
        cout<<"feeId "<<feeId<<endl;
        cout<<"chId "<<chId<<endl;
        cout<<"adcRange "<<adcRange<<endl;
        cout<<"infoCode "<<infoCode<<endl;
    }
    void Clear(){
        timestampLsb=0xFFFFFFFF;
        infoField=0xFFFFFFFF;
        adcData=0xFFFFFFFF;
        dataType=0xFFFFFFFF;
        feeId=0xFFFFFFFF;
        chId=0xFFFFFFFF;
        adcRange=0xFFFFFFFF;
        infoCode=0xFFFFFFFF;
        sampleLengh=0xFFFFFFFF;
    }

};

class rawaida_info {
public :
   // Declaration of leaf types
   Long64_t        evt;
   Long64_t        timestamp;
   Long64_t        extTimestamp;
   Int_t           feeNo;
   Int_t           chNo;
   Int_t           dssdNo;
   Int_t           stripNo;
   Short_t         infoCode;
   Short_t         rangeType;
   Int_t           adcData;

   rawaida_info(){}
   virtual ~rawaida_info(){}
   void     Show(){
       cout<<"timestamp "<<timestamp<<endl;
       cout<<"external timestamp "<<extTimestamp<<endl;
       cout<<"info code "<<infoCode<<endl;
       cout<<"fee no. "<<feeNo<<endl;
       cout<<"channel no. "<<chNo<<endl;
       cout<<"strip no. "<<stripNo<<endl;
       cout<<"dssd no. "<<dssdNo<<endl;
       cout<<"range type "<<rangeType<<endl;
   }

   void  Clear(){
       adcData=-9999;
       timestamp=-9999;
       extTimestamp=-9999;
       infoCode=-9999;
       feeNo=-9999;
       chNo=-9999;
       stripNo=-9999;
       dssdNo=-9999;
       rangeType=-9999;
   }
};

#endif


