#ifndef __BEAM_H
#define __BEAM_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"

#define kMaxBeamInfo 1
#define kMaxRIPSInfo 1
#define kMaxTOFInfo 1



//using namespace std;
/*!
  Container for the full beam, tof, beta and pid information
*/
class Beam : public TObject {
public:
  //! default constructor
  Beam(){
    Clear();
  }
  virtual ~Beam(){}
  //! Clear all information
  void Clear(Option_t *option = ""){
    fts = 0;
    memset(faoq,0.,sizeof(faoq));
    memset(faoqc,0.,sizeof(faoqc));
    memset(fzet,0.,sizeof(fzet));
    memset(fzetc,0.,sizeof(fzetc));
    memset(ftof,0.,sizeof(ftof));
    memset(fbeta,0.,sizeof(fbeta));
    memset(fdelta,0.,sizeof(fdelta));
    /*
    for(unsigned short j=0;j<6;j++){
      faoq[j] = sqrt(-1.);
      faoqc[j] = sqrt(-1.);
      fzet[j] = sqrt(-1.);
      fzetc[j] = sqrt(-1.);
    }

    for(unsigned short j=0;j<3;j++){
      ftof[j] = sqrt(-1.);
      fbeta[j] = sqrt(-1.);
    }
    for(unsigned short j=0;j<4;j++){
      fdelta[j] = sqrt(-1.);
    }
    */
  }
  virtual void Copy(Beam& obj){
      obj.SetTimestamp(fts);
      for (Int_t i=0;i<kMaxBeamInfo;i++){
          obj.SetAQ(i,faoq[i]);
          obj.SetCorrAQ(i,faoqc[i]);
          obj.SetZ(i,fzet[i]);
          obj.SetCorrZ(i,fzetc[i]);
          obj.SetBeta(i,fbeta[i]);
      }
      for (Int_t i=0;i<kMaxTOFInfo;i++){
          obj.SetTOF(i,ftof[i]);
      }
      for (Int_t i=0;i<kMaxRIPSInfo;i++){
          obj.SetDelta(i,fdelta[i]);
      }
  }

  //! Set the timestamp
  void SetTimestamp(unsigned long long ts){fts = ts;}

  //! Set the A/Q ratio
  void SetAQ(unsigned short j, double aoq){
    if( j>kMaxBeamInfo-1) return;
    faoq[j] = aoq;
  }
  //! Set the A/Q corrected ratio
  void SetCorrAQ(unsigned short j, double aoq){
    if( j>kMaxBeamInfo-1) return;
    faoqc[j] = aoq;
  }
  //! Set the Z number
  void SetZ(unsigned short j, double zet){
    if( j>kMaxBeamInfo-1) return;
    fzet[j] = zet;
  }
  //! Set the Z corrected
  void SetCorrZ(unsigned short j, double zet){
    if( j>kMaxBeamInfo-1) return;
    fzetc[j] = zet;
  }
  //! Set both A/Q and Z
  void SetAQZ(unsigned short j, double aoq, double zet){
    if( j>kMaxBeamInfo-1) return;
    faoq[j] = aoq;
    fzet[j] = zet;
  }
  //! Set the time-of-flight
  void SetTOF(unsigned short j, double tof){
    if( j>kMaxTOFInfo-1) return;
    ftof[j] = tof;
  }
  //! Set the beta
  void SetBeta(unsigned short j, double beta){
    if( j>kMaxBeamInfo-1) return;
    fbeta[j] = beta;
  }

  //! Set the delta
  void SetDelta(unsigned short j, double delta){
    if( j>kMaxRIPSInfo-1) return;
    fdelta[j] = delta;
  }

  //! Correct the A/Q ratio based on position
  void CorrectAQ(unsigned short j, double corr){
    if( j>kMaxBeamInfo-1) return;
    faoqc[j] = faoq[j] + corr;
  }

  //! Get timestamp
  unsigned long long GetTimestamp(){return fts;}

  //! Get the A/Q ratio
  double GetAQ(unsigned short j){
    if( j>kMaxBeamInfo-1) return sqrt(-1.);
    return faoq[j];
  }
  //! Get the corrected A/Q ratio
  double GetCorrAQ(unsigned short j){
    if( j>kMaxBeamInfo-1) return sqrt(-1.);
    return faoqc[j];
  }
  //! Get the Z number
  double GetZ(unsigned short j){
    if( j>kMaxBeamInfo-1) return sqrt(-1.);
    return fzet[j];
  }
  //! Get the time-of-flight
  double GetTOF(unsigned short j){
    if( j>2) return sqrt(-1.);
    return ftof[j];
  }
  //! Get beta
  double GetBeta(unsigned short j){
    if( j>2) return sqrt(-1.);
    return fbeta[j];
  }
  //! Get Delta
  double GetDelta(unsigned short j){
    if( j>2) return sqrt(-1.);
    return fdelta[j];
  }

protected:
  unsigned long long fts;
  //! A/Q for 3-5, 5-7, 3-7
  double faoq[kMaxBeamInfo];
  //! corrected A/Q
  double faoqc[kMaxBeamInfo];
  //! Z for 3-5, 5-7, 3-7
  double fzet[kMaxBeamInfo];
  //! corrected Z
  double fzetc[kMaxBeamInfo];

  //! time-of-flight for 3-5, 5-7, 3-7
  double ftof[kMaxTOFInfo];
  //! beta for 3-7, 8-11, 7-8
  double fbeta[kMaxBeamInfo];
  //! delta momentum 3-5, 5-7, 3-7
  double fdelta[kMaxRIPSInfo];

  /// \cond CLASSIMP
  ClassDef(Beam,1);
  /// \endcond
};


#endif
