#ifndef AIDA_H
#define AIDA_H
#include <iostream>
#include <vector>
#include <map>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "AIDAdefs.h"
using namespace std;

class AIDASimpleStruct : public TObject {
public:
  //! default constructor
  AIDASimpleStruct(){
    Clear();
  }
  virtual ~AIDASimpleStruct(){}

  //! Clear the music information
  virtual void Clear(){
      fid=0.;
      fts=0;
      fevt=0;
      fx=-9999.;
      fy=-9999.;
      fz=-1;
      fex=-9999.;
      fey=-9999.;
      fncx=0;
      fncy=0;
      fnx=0;
      fny=0;
      fnz=0;
      fniz=0;
      frflag=0;
      fdrflag=0;
      ftw=-9999.;
      fdtion=-9999.;
      fdz=0;
      fmult=0;
      fmindy=9999;
      fmindx=9999;
      fsumexyrank=0;
      fmaxelastdssd=-9999;
      fxts=0;
      fyts=0;
      fxtsmax=0;
      fytsmax=0;
      fminx = -11;
      fmaxx = -11;
      fminy = -11;
      fmaxy = -11;
      ftdiffmin= -99999;
      ftdiffmax= -99999;
  }

  //! copy cluster
  virtual void Copy(AIDASimpleStruct& obj){
      obj.SetID(fid);
      obj.SetEventNumber(fevt);
      obj.SetTimestamp(fts);
      obj.SetMult(fmult);
      obj.SetHitPosition(fx,fy,fz);
      obj.SetXEnergy(fex);
      obj.SetYEnergy(fey);
      obj.SetXClusterMult(fncx);
      obj.SetYClusterMult(fncy);
      obj.SetXMult(fnx);
      obj.SetYMult(fny);
      obj.SetZMult(fnz);
      obj.SetStripMultFlag(fniz);
      obj.SetRankingFlag(frflag);
      obj.SetEDiffRankingFlag(fdrflag);
      obj.SetTimeWidth(ftw);
      obj.SetDtIon(fdtion);
      obj.SetDZ(fdz);
      obj.SetMinimumDistanceX(fmindx);
      obj.SetMinimumDistanceY(fmindy);
      obj.SetSumEXYRank(fsumexyrank);
      obj.SetMaxELastDSSD(fmaxelastdssd);
      obj.SetXTimestamp(fxts);
      obj.SetYTimestamp(fyts);
      obj.SetXMaxTimestamp(fxtsmax);
      obj.SetYMaxTimestamp(fytsmax);
      obj.SetXMinPos(fminx);
      obj.SetXMaxPos(fmaxx);
      obj.SetYMinPos(fminy);
      obj.SetYMaxPos(fmaxy);
      obj.SetTimeDifferenceMax(ftdiffmax);
      obj.SetTimeDifferenceMin(ftdiffmin);
  }

  //! Set XY
  void SetHitPosition(double x,double y,short z){
      fx=x;
      fy=y;
      fz=z;
  }

  //! xy area of cluster
  void SetXMinPos(short minx) {fminx=minx;}
  void SetXMaxPos(short maxx) {fmaxx=maxx;}
  void SetYMinPos(short miny) {fminy=miny;}
  void SetYMaxPos(short maxy) {fmaxy=maxy;}

  //! Set ID of evnt
  void SetID(unsigned char id){fid = id;}

  //! Set Evt number
  void SetEventNumber(unsigned int evt){fevt = evt;}

  //! Set the X strips sum energy
  void SetXEnergy(double xenergy){fex = xenergy;}
  //! Set the Y strips sum energy
  void SetYEnergy(double yenergy){fey = yenergy;}

  //! Set the X strips ealiest timestamp
  void SetXTimestamp(unsigned long long xts){fxts = xts;}
  //! Set the Y strips ealiest timestamp
  void SetYTimestamp(unsigned long long yts){fyts = yts;}

  //! Set the X strips latest time
  void SetXMaxTimestamp(double xts){fxtsmax = xts;}


  //! Set time diffrernce between X-Y
  void SetTimeDifferenceMin(int tdiffmin){ftdiffmin = tdiffmin;}
  void SetTimeDifferenceMax(int tdiffmax){ftdiffmax = tdiffmax;}


  //! Set the Y strips latest time
  void SetYMaxTimestamp(double yts){fytsmax = yts;}

  //! Set the event multiplicity
  void SetMult(unsigned short mult){fmult = mult;}

  //! Set the X strips multiplicity
  void SetXMult(unsigned short xmult){fnx = xmult;}
  //! Set the Y strips multiplicity
  void SetYMult(unsigned short ymult){fny = ymult;}

  //! Set the X strips multiplicity
  void SetXClusterMult(unsigned short xcmult){fncx = xcmult;}
  //! Set the Y strips multiplicity
  void SetYClusterMult(unsigned short ycmult){fncy = ycmult;}


  //! Set the DSSD  multiplicity
  void SetZMult(unsigned short zmult){fnz = zmult;}

  //! Set the timestamp
  void SetTimestamp(unsigned long long ts){fts = ts;}

  //! Set Energy ratio ranking flag
  void SetRankingFlag(unsigned char rankingflag){frflag = rankingflag;}

  //! Set Energy difference ranking flag
  void SetEDiffRankingFlag(unsigned char rankingflag){fdrflag = rankingflag;}

  //! Set Minimum distance X
  void SetMinimumDistanceX(double mindx){fmindx=mindx;}
  //! Set Minimum distance Y
  void SetMinimumDistanceY(double mindy){fmindy=mindy;}

  //! Set Multiple hit in strip flag
  void SetStripMultFlag(unsigned short niz){fniz = niz;}

  //! Set Event Time Width
  void SetTimeWidth(double tw){ftw=tw;}

  //! Set time distance to ion
  void SetDtIon(double dtion){fdtion=dtion;}

  //! Set Z correction
  void SetDZ(unsigned short dz){fdz=dz;}

  //! Set EX+EY rank
  void SetSumEXYRank(unsigned short sumexyrank){fsumexyrank = sumexyrank;}

  //! Set Max Energy of last layer
  void SetMaxELastDSSD(double maxelastdssd){fmaxelastdssd = maxelastdssd;}

  //! Get ID of evnt
  unsigned char GetID(){return fid;}

  unsigned int GetEventNumber(){return fevt;}

  double GetHitPositionX(){return fx;}
  double GetHitPositionY(){return fy;}
  short GetHitPositionZ(){return fz;}

  //! Get xy area of cluster
  short GetMinHitPositionX(){return fminx;}
  short GetMaxHitPositionX(){return fmaxx;}
  short GetMinHitPositionY(){return fminy;}
  short GetMaxHitPositionY(){return fmaxy;}

  //! Get the X strips sum energy
  double GetXEnergy(){return fex;}
  //! Get the Y strips sum energy
  double GetYEnergy(){return fey;}

  //! Get the X strips ealiest timestamp
  unsigned long long GetXMinTimestamp(){return fxts;}
  //! Get the Y strips ealiest timestamp
  unsigned long long GetYMinTimestamp(){return fyts;}

  //! Get the X strips latest timestamp
  unsigned long long GetXMaxTimestamp(){return fxtsmax;}
  //! Get the Y strips latest timestamp
  unsigned long long GetYMaxTimestamp(){return fytsmax;}

  //! Get the time diffrerence of Y and X strips (min)
  int GetTimeDifferenceMin(){return ftdiffmin;}
  //! Get the time diffrerence of Y and X strips (max)
  int GetTimeDifferenceMax(){return ftdiffmax;}


  //! Get the event multiplicity
  unsigned short GetMultiplicity(){return fmult;}


  //! Get the X clusters multiplicity
  unsigned short GetXClusterMultiplicity(){return fncx;}
  //! Get the Y clusters multiplicity
  unsigned short GetYClusterMultiplicity(){return fncy;}

  //! Get the X strips multiplicity
  unsigned short GetXMultiplicity(){return fnx;}
  //! Get the Y strips multiplicity
  unsigned short GetYMultiplicity(){return fny;}

  //! Get the Z strips multiplicity
  unsigned short GetZMultiplicity(){return fnz;}

  //! Get the timestamp
  unsigned long long GetTimestamp(){return fts;}

  //! Get Energy ratio ranking flag
  unsigned char GetRankingFlag(){return frflag;}

  //! Get Energy difference ranking flag
  unsigned char GetEDiffRankingFlag(){return fdrflag;}

  //! Get minimum cluster distance X
  double GetMinimumDistanceX(){return fmindx;}
  //! Get minimum cluster distance Y
  double GetMinimumDistanceY(){return fmindy;}

  //! Get Multiple hit in strip flag
  unsigned short GetStripMultFlag(){return fniz;}

  //! Get time distance to ion
  double GetDtIon(){return fdtion;}

  //! Get even time width
  double GetTimeWidth(){return ftw;}

  //! Get Multiple hit in strip flag
  unsigned short GetDZ(){return fdz;}

  //! Get ex+ey rank
  unsigned short GetSumEXYRank(){return fsumexyrank;}

  //! Set Max Energy of last layer
  double GetMaxELastDSSD(){return fmaxelastdssd;}


  //! Printing information
  void Print(Option_t *option = "") const {
    cout <<"\tts "<< fts;
    cout << "\tX " << fx;
    cout << "\tY " << fy;
    cout << "\tZ " << fz;
    cout << "\tX energy " << fex;
    cout << "\tY energy " << fey;
    cout << "\tX multiplicity " << fnx;
    cout << "\tY multiplicity " << fny;
    cout << "\ttimestamp " << fts;
    cout << "\tranking flag " << frflag;
    return;
  }

protected:
  unsigned long long fts;

  unsigned int fevt;

  //! id :ion or beta
  unsigned char fid;

  //! event multiplicity
  unsigned short fmult;

  //! DSSDs coordinator
  double fx;
  double fy;
  short fz;

  //! area of cluster
  short fminx;
  short fmaxx;
  short fminy;
  short fmaxy;

  //! store the hit number
  //vector<short*> fhitsno;

  //! the energy lab system
  double fex;
  double fey;
  //! xy hit multiplicity
  unsigned short fnx;
  unsigned short fny;
  //! xy cluster multiplicity
  unsigned short fncx;
  unsigned short fncy;
  //! number of clusters
  unsigned short fnz;
  //! number of event with mult>1 for 1 dssd
  unsigned short fniz;
  //! ranking flag
  unsigned char frflag;
  //! ranking flag EY-EX
  unsigned char fdrflag;

  //!ealiest time stamp in x
  unsigned long long fxts;
  //!ealiest time stamp in y
  unsigned long long fyts;
  //!lastest time stamp in x
  unsigned long long fxtsmax;
  //!lastest time stamp in y
  unsigned long long fytsmax;

  //! time diffrence between y and x cluster
  int ftdiffmin;
  int ftdiffmax;

  //! ex+ey rank
  unsigned short fsumexyrank=0;

  //! minimum cluster distance x
  double fmindx;
  //! minimum cluster distance y
  double fmindy;

  //! event time width in us
  double ftw;
  //! time distance with prev ion in us
  double fdtion;

  //! z correction (if applicable)
  unsigned short fdz;

  //! max energy of the last layer (for light particle veto method)
  double fmaxelastdssd;

  /// \cond CLASSIMP
  ClassDef(AIDASimpleStruct,1);
  /// \endcond
};



class AIDAHit : public TObject {
public:
  //! default constructor
  AIDAHit(){
    Clear();
  }
  virtual ~AIDAHit(){}
  //! constructor with individual values
  AIDAHit(short range, short fee, short ch, short id, int xy, int z, double en, int adc,unsigned short hitsadded, unsigned long long int ts){
    ffee = fee;
    fch = ch;
    fid = id;
    fxy=xy;
    fz=z;
    fadc=adc;
    fen = en;
    fhitsadded = hitsadded;
    fts = ts;
    frange = range;
  }
  //! Clear the music information
  virtual void Clear(Option_t *option = ""){
      fch = -1;
      ffee = -1;
      fid = -1;
      fxy = -1;
      fz = -1;
      fadc = -9999;
      fen = -9999.;
      fhitsadded = 0;
      fts = 0;
      ffastts = 0;
      frange = 0;
  }

  //! Copy hits
  virtual void Copy(AIDAHit& obj){
      obj.SetFEE(ffee);
      obj.SetFEEChannel(fch);
      obj.SetID(fid);
      obj.SetXY(fxy);
      obj.SetZ(fz);
      obj.SetADC(fadc);
      obj.SetEnergy(fen);
      obj.SetTimestamp(fts);
      obj.SetFastTimestamp(ffastts);
      obj.SetRange(frange);
  }

  //! Set energy range (low energy 0 or high energy 1)
  void SetRange(short range){frange = range;}

  //! Set the strip ID
  void SetID(short id){fid = id;}
  //! Set XY
  void SetXY(short xy){fxy = xy;}
  //! Set Z
  void SetZ(short z){fz= z;}
  //! Set hit position in DSSD
  void SetStrip(short xy,short z){
      fxy = xy;
      fz = z;
  }

  //! Set FEE number
  void SetFEE(short fee){ffee = fee;}
  //! Set channel number
  void SetFEEChannel(short ch){fch = ch;}

  //! Set the energy
  void SetEnergy(double energy){fen = energy;}

  //! Set the raw ADC value
  void SetADC(int adc){fadc = adc;}
  //! Set the timestamp
  void SetTimestamp(unsigned long long int ts){fts = ts;}
  //! Set the fast timestamp
  void SetFastTimestamp(unsigned long long int fts){ffastts = fts;}

  //! Set current hits
  void SetHitsAdded(unsigned short hitsadded){fhitsadded = hitsadded;}

  //! Get the range
  short GetRange(){return frange;}
  //! Get the ID
  short GetID(){return fid;}
  //! Get XY position in DSSD
  short GetXY(){return fxy;}
  //! Get Z position in DSSD
  short GetZ(){return fz;}

  //! Get Fee number
  short GetFEE(){return ffee;}
  //! Get Fee channel number
  short GetFEEChannel(){return fch;}

  //! Get the energy
  double GetEnergy(){return fen;}

  //! Get the timestamp
  unsigned long long int GetTimestamp(){return fts;}
  //! Get the fast timestamp
  unsigned long long int GetFastTimestamp(){return ffastts;}

  //! Get the raw ADC value
  int GetADC(){return fadc;}
  //! Get current hits
  unsigned short GetHitsAdded(){return fhitsadded;}

  //! Printing information
  void Print(Option_t *option = "") const {
    cout << "ID " << fid;
    cout << "\tDSSD " << fz;
    cout << "\tStrip " << fxy;
    cout << "\tadc " << fadc;
    cout << "\tenergy " << fen;
    cout << "\ttimestamp " << fts;
    cout << "\thits added " << fhitsadded << endl;
    return;
  }

protected:
  //! energy range : (0: high gain, 1: low gain)
  short frange;
  //! translate into channel ID number
  short fid;
  //! translate into the DSSDs coordinator
  short fxy;
  short fz;
  //! FEE information
  short ffee;
  short fch;

  //! the energy lab system
  double fen;
  //! the raw adc value
  int fadc;
  //! the timestamp
  unsigned long long fts;
  //! the fast timestamp
  unsigned long long ffastts;

  //! current hits
  unsigned short fhitsadded;

  /// \cond CLASSIMP
  ClassDef(AIDAHit,1);
  /// \endcond
};


class AIDACluster : public TObject {
public:
  //! default constructor
  AIDACluster(){
    Clear();
  }
  virtual ~AIDACluster(){}
  //! constructor with individual values
  AIDACluster(unsigned short x, unsigned short y, unsigned short z, unsigned short multx,unsigned short multy,
              double xenergy,double yenergy,unsigned short clustersadded, unsigned long long int ts, unsigned char rankingflag,unsigned char ediffrankingflag, unsigned short sumexyrank)
  {
    fpos.SetXYZ(x,y,z);
    fsumenx=xenergy;
    fsumeny=yenergy;
    fnx=multx;
    fny=multy;
    fclustersadded = clustersadded;
    fcts = ts;    
    frankingflag = rankingflag;
    fediffrankingflag = ediffrankingflag;
    fsumexyrank = sumexyrank;
  }
  //! Clear the music information
  virtual void Clear(){
      fpos.SetXYZ(-1,-1,-1);
      fposMAXE.SetXYZ(-1,-1,-1);
      fsumenx=-9999.;
      fsumeny=-9999.;
      fcts=0;
      fclustersadded=0;
      fnx=0;
      fny=0;
      fcfastts=0;
      frankingflag = 0;
      fediffrankingflag=0;
      fsumexyrank=0;
      fxts = 0;
      fyts = 0;
      fxtsmax = 0;
      fytsmax = 0;
      fminx = -11;
      fmaxx = -11;
      fminy = -11;
      fmaxy = -11;
      ftdiffmin= -99999;
      ftdiffmax= -99999;
  }

  //! copy cluster
  virtual void Copy(AIDACluster& obj){
      obj.SetTimestamp(fcts);
      obj.SetHitPosition(fpos.X(),fpos.Y(),fpos.Z());
      obj.SetHitPositionMaxE(fposMAXE.X(),fposMAXE.Y(),fposMAXE.Z());
      obj.SetXEnergy(fsumenx);
      obj.SetYEnergy(fsumeny);
      obj.SetXMult(fnx);
      obj.SetYMult(fny);
      obj.SetFastTimestamp(fcfastts);
      obj.SetClustersAdded(fclustersadded);
      obj.SetRankingFlag(frankingflag);
      obj.SetEDiffRankingFlag(fediffrankingflag);
      obj.SetSumEXYRank(fsumexyrank);
      obj.SetXTimestamp(fxts);
      obj.SetYTimestamp(fyts);
      obj.SetXMaxTimestamp(fxtsmax);
      obj.SetYMaxTimestamp(fytsmax);
      obj.SetXMinPos(fminx);
      obj.SetXMaxPos(fmaxx);
      obj.SetYMinPos(fminy);
      obj.SetYMaxPos(fmaxy);
      obj.SetTimeDifferenceMax(ftdiffmax);
      obj.SetTimeDifferenceMin(ftdiffmin);
  }

  //! Set XY
  void SetHitPosition(double x,double y,double z){fpos.SetXYZ(x,y,z);}
  void SetHitPositionMaxE(double x,double y,double z){fposMAXE.SetXYZ(x,y,z);}

  //! xy area of cluster
  void SetXMinPos(short minx) {fminx=minx;}
  void SetXMaxPos(short maxx) {fmaxx=maxx;}
  void SetYMinPos(short miny) {fminy=miny;}
  void SetYMaxPos(short maxy) {fmaxy=maxy;}

  //! Set the X strips sum energy
  void SetXEnergy(double xenergy){fsumenx = xenergy;}
  //! Set the Y strips sum energy
  void SetYEnergy(double yenergy){fsumeny = yenergy;}

  //! Set the X strips earliest time
  void SetXTimestamp(double xts){fxts = xts;}

  //! Set the Y strips earliest time
  void SetYTimestamp(double yts){fyts = yts;}   

  //! Set the X strips latest time
  void SetXMaxTimestamp(double xts){fxtsmax = xts;}

  //! Set the Y strips latest time
  void SetYMaxTimestamp(double yts){fytsmax = yts;}

  //! Set time diffrernce between X-Y
  void SetTimeDifferenceMin(int tdiffmin){ftdiffmin = tdiffmin;}
  void SetTimeDifferenceMax(int tdiffmax){ftdiffmax = tdiffmax;}

  void CalculateTimeDifference(){
      //! tdiff x min
      int tdiffyminxmin = abs((long long)fyts-(long long)fxts);
      int tdiffymaxxmin = abs((long long)fytsmax-(long long)fxts);
      if (tdiffyminxmin > tdiffymaxxmin){
          ftdiffmin = tdiffymaxxmin;
      }else{
          ftdiffmin = tdiffyminxmin;
      }
      //! tdiff x max
      int tdiffyminxmax = abs((long long)fyts-(long long)fxtsmax);
      int tdiffymaxxmax = abs((long long)fytsmax-(long long)fxtsmax);
      if (tdiffyminxmax > tdiffymaxxmax){
          ftdiffmax = tdiffymaxxmax;
      }else{
          ftdiffmax = tdiffyminxmax;
      }
  }


  //! Set the X strips multiplicity
  void SetXMult(unsigned short xmult){fnx = xmult;}
  //! Set the Y strips multiplicity
  void SetYMult(unsigned short ymult){fny = ymult;}


  //! Set the timestamp
  void SetTimestamp(unsigned long long int ts){fcts = ts;}

  //! Set the fast timestamp
  void SetFastTimestamp(unsigned long long int fts){fcfastts = fts;}

  //! Set current cluster number
  void SetClustersAdded(unsigned short clustersadded){fclustersadded = clustersadded;}

  //! Set XY energy diffrence
  void SetRankingFlag(unsigned char rankingflag){frankingflag = rankingflag;}

  //! Set Energy difference ranking flag
  void SetEDiffRankingFlag(unsigned char rankingflag){fediffrankingflag = rankingflag;}

  //! Set Sum energy of X and Y rank
  void SetSumEXYRank(unsigned short sumexyrank){fsumexyrank = sumexyrank;}


  //! Get Hit position
  TVector3 GetHitPosition(){return fpos;}
  Double_t GetHitPositionX(){return fpos.X();}
  Double_t GetHitPositionY(){return fpos.Y();}
  Double_t GetHitPositionZ(){return fpos.Z();}


  //! Get Hit max energy
  TVector3 GetHitPositionMaxE(){return fposMAXE;}
  Double_t GetHitPositionMaxEX(){return fposMAXE.X();}
  Double_t GetHitPositionMaxEY(){return fposMAXE.Y();}
  Double_t GetHitPositionMaxEZ(){return fposMAXE.Z();}


  //! Get xy area of cluster
  short GetMinHitPositionX(){return fminx;}
  short GetMaxHitPositionX(){return fmaxx;}
  short GetMinHitPositionY(){return fminy;}
  short GetMaxHitPositionY(){return fmaxy;}


  //! Get the X strips sum energy
  double GetXEnergy(){return fsumenx;}
  //! Get the Y strips sum energy
  double GetYEnergy(){return fsumeny;}

  //! Get the X strips ealiest timestamp
  unsigned long long GetXTimestamp(){return fxts;}
  //! Get the Y strips ealiest timestamp
  unsigned long long GetYTimestamp(){return fyts;}

  //! Get the X strips latest timestamp
  unsigned long long GetXMaxTimestamp(){return fxtsmax;}
  //! Get the Y strips latest timestamp
  unsigned long long GetYMaxTimestamp(){return fytsmax;}

  //! Get the time diffrerence of Y and X strips (min)
  int GetTimeDifferenceMin(){return ftdiffmin;}
  //! Get the time diffrerence of Y and X strips (max)
  int GetTimeDifferenceMax(){return ftdiffmax;}



  //! Get the X strips multiplicity
  unsigned short GetXMultiplicity(){return fnx;}
  //! Get the Y strips multiplicity
  unsigned short GetYMultiplicity(){return fny;}

  //! Get the timestamp
  unsigned long long int GetTimestamp(){return fcts;}

  //! Get the fast timestamp
  unsigned long long int GetFastTimestamp(){return fcfastts;}

  //! Get current cluster number
  unsigned short GetClustersAdded(){return fclustersadded;}

  //! Get Energy diffrence ranking flag
  unsigned char GetRankingFlag(){return frankingflag;}

  //! Get Energy difference ranking flag
  unsigned char GetEDiffRankingFlag(){return fediffrankingflag;}

  //! Get Sum energy ranking
  unsigned short GetSumEXYRank(){return fsumexyrank;}

  //! Printing information
  void Print(Option_t *option = "") const {
    cout << "Cluster No. " << fclustersadded;
    cout << "\tX " << fpos.X();
    cout << "\tY " << fpos.Y();
    cout << "\tZ " << fpos.Z();
    cout << "\tX energy " << fsumenx;
    cout << "\tY energy " << fsumeny;
    cout << "\tX multiplicity " << fnx;
    cout << "\tY multiplicity " << fny;
    cout << "\ttimestamp " << fcts;
    cout << "\tranking flag " << frankingflag;
    return;
  }

protected:
  //! translate into the DSSDs coordinator
  TVector3 fpos;

  TVector3 fposMAXE;

  //! store the hit number
  //vector<short*> fhitsno;

  //! the energy lab system
  double fsumenx;
  double fsumeny;

  //! number of hit in one cluster
  unsigned short fnx;
  unsigned short fny;

  //! x - y time stamp
  unsigned long long fxts;
  unsigned long long fyts;
  unsigned long long fxtsmax;
  unsigned long long fytsmax;

  //! time diffrence between y and x cluster
  int ftdiffmin;
  int ftdiffmax;

  //! area of cluster
  short fminx;
  short fmaxx;
  short fminy;
  short fmaxy;

  //! the ealiest timestamp
  unsigned long long fcts;

  //! the ealiest fast timestamp
  unsigned long long fcfastts;

  //! current hits
  unsigned short fclustersadded;

  //! ranking mode: if cluster satisfy ranking mode condition = 1, otherwise =0
  unsigned char frankingflag;

  unsigned char fediffrankingflag;


  unsigned short fsumexyrank;

  /// \cond CLASSIMP
  ClassDef(AIDACluster,1);
  /// \endcond
};

class AIDA : public TObject
{
public:
    //! default constructor
    AIDA(){

      //for (Int_t i=0;i<NumDSSD;i++){
      //    for (Int_t j=0;j<NumStrXY;j++){
      //        fdssd_thr[i][j] = 0.;
      //        fdssd_cal[i][j][0] = 0.;
      //        fdssd_cal[i][j][1] = 1.;
      //    }
      //}
      Clear();
    }
    virtual ~AIDA(){}
    //! Clear the AIDA information
    virtual void Clear(Option_t *option = ""){
      faidats = 0;
      faidatw = 0;
      ftsprevion = 0;
      fmult = 0;
      fnclusters = 0;
      fmaxz = 0;
      fclustermultz = 0;
      fhitmultz = 0;
      fdmaxz = 0;
      /*
      for (int i=0;i<NumDSSD;i++){
          fmultx[i]=0;
          fmulty[i]=0;
      }
      */
      memset(fmultx,0,sizeof(fmultx));
      memset(fmulty,0,sizeof(fmulty));

      memset(fsumx,0,sizeof(fsumx));
      memset(fsumy,0,sizeof(fsumy));

      memset(fnhitsz,0,sizeof(fnhitsz));

      memset(fnclustersz,0,sizeof(fnclustersz));

      //! 2017 july 26 added
      memset(fstripsmult,0,sizeof(fstripsmult));
      memset(fstripsmultx1,0,sizeof(fstripsmultx1));
      memset(fstripsmultx2,0,sizeof(fstripsmultx2));
      memset(fstripsmulty1,0,sizeof(fstripsmulty1));
      memset(fstripsmulty2,0,sizeof(fstripsmulty2));

      memset(fnxclustersz,0,sizeof(fnxclustersz));
      memset(fnyclustersz,0,sizeof(fnyclustersz));

      memset(fmindx,0,sizeof(fmindx));
      memset(fmindy,0,sizeof(fmindy));


      //! Dealocating memory
      for (size_t idx=0;idx<fhits.size();idx++){
          delete fhits[idx];
      }
      fhits.clear();

      for (size_t idx=0;idx<fclusters.size();idx++){
          delete fclusters[idx];
      }
      fclusters.clear();

      ftype=0;

    }

    void Copy(AIDA& obj){
        for(vector<AIDACluster*>::iterator cluster=fclusters.begin(); cluster!=fclusters.end(); cluster++){
          AIDACluster* clonecluster = new AIDACluster;
          AIDACluster* origincluster = *cluster;
          origincluster->Copy(*clonecluster);
          //! CALIBRATE TIME HERE!
          //clonecluster->SetTimestamp(clonecluster->GetTimestamp()*ClockResolution);

      //if (!(clonecluster->GetHitPositionZ()==1&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()<64))
      //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>125))
      //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>50&&clonecluster->GetHitPositionY()<75))
      //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>98&&clonecluster->GetHitPositionY()<111))
      //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()>126))
          obj.AddCluster(clonecluster);
        }
        //obj.SetTimestamp(faidats*ClockResolution);
        obj.SetTimestamp(faidats);
        obj.SetType(ftype);
        obj.SetMult(fmult);
        obj.SetNHitZ(fnhitsz);
        obj.SetHitMultZ(fhitmultz);
        obj.SetMultX(fmultx);
        obj.SetMultY(fmulty);
        obj.SetSumX(fsumx);
        obj.SetSumY(fsumy);
        //obj.SetNClusters(fnclusters);
        //obj.SetNClusterZ(fnclustersz);
        obj.SetMaxZ(fmaxz);
        //obj.SetClusterMultZ(fclustermultz);
    }

    void SetTimestamp(unsigned long long ts){faidats = ts;}
    //! added 26th July 2017
    void SetTimeWindow(unsigned long long tw){faidatw = tw;}

    void  SetPrevIonTimestamp(unsigned long long ts){ftsprevion = ts;}

    //!clear all hits
    void ClearAllHits(){
        //! Dealocating memory also
        for (size_t idx=0;idx<fhits.size();idx++){
            delete fhits[idx];
        }
        fhits.clear();
    }

    void ClearAllClusters(){
        //! Dealocating memory also
        for (size_t idx=0;idx<fclusters.size();idx++){
            delete fclusters[idx];
        }
        fclusters.clear();
    }

    //! Add a hit
    void AddHit(AIDAHit* hit){
      //!newly added
      hit->SetHitsAdded(fmult);
      fhits.push_back(hit);
      //! newly added
      Int_t z=(Int_t)hit->GetZ();
      if (hit->GetXY() < NumStrX) {
          fmultx[z]++;
          fsumx[z]+=hit->GetEnergy();
      }else {
          fmulty[z]++;
          fsumy[z]+=hit->GetEnergy();
      }

      int Z=(int) hit->GetZ();
      fnhitsz[Z]++;
      if (fnhitsz[Z]==1) fhitmultz++;
      fmult++;
    }

    //! Add more hits
    void AddHits(vector<AIDAHit*> hits){
      fmult += hits.size();
      for(vector<AIDAHit*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
          //set hit add here!
          //
          fhits.push_back(*hit);
      }
    }

    //! Set all clusters
    void SetHits(vector<AIDAHit*> hits){
      fmult = hits.size();
      fhits = hits;
    }


    //! Set multiplicity
    void SetMult(unsigned short mult){
        fmult =  mult;
    }

    //! set number of dssd with at least 1 hit!
    void SetHitMultZ(unsigned short hitmultz){
        fhitmultz = hitmultz;
    }

    //! Set the nhits in z
    void SetNHitZ(unsigned short nhitsz[NumDSSD]){
        memcpy(fnhitsz,nhitsz,NumDSSD*sizeof(unsigned short));
    }

    //! Set the X strip multiplicity of the event
    void SetMultX(unsigned short multx[NumDSSD]){
        memcpy(fmultx,multx,NumDSSD*sizeof(unsigned short));
    }

    //! Set the Y strip multiplicity of the event
    void SetMultY(unsigned short multy[NumDSSD]){
        memcpy(fmulty,multy,NumDSSD*sizeof(unsigned short));
    }

    //! Set the X strip sum energy of the event
    void SetSumX(double sumx[NumDSSD]){
        memcpy(fsumx,sumx,NumDSSD*sizeof(double));
    }

    //! Set the Y strip sum energy of the event
    void SetSumY(double sumy[NumDSSD]){
        memcpy(fsumy,sumy,NumDSSD*sizeof(double));
    }

    //! Set per silicon channel multi hit flag
    void AddStripMult(unsigned short dssd){
        fstripsmult[dssd]++;
    }

    void SetStripMultX1Flag(unsigned short dssd,unsigned long long val){
        fstripsmultx1[dssd]=val;
    }
    void SetStripMultX2Flag(unsigned short dssd,unsigned long long val){
        fstripsmultx2[dssd]=val;
    }
    void SetStripMultY1Flag(unsigned short dssd,unsigned long long val){
        fstripsmulty1[dssd]=val;
    }
    void SetStripMultY2Flag(unsigned short dssd,unsigned long long val){
        fstripsmulty2[dssd]=val;
    }

    //! Set per silicon channel multi hit flag
    void SetStripMultX1FlagMask(unsigned short dssd,unsigned short stripNo){
        fstripsmultx1[dssd]=fstripsmultx1[dssd]|((ULong64_t)0x1<<stripNo);
    }
    void SetStripMultX2FlagMask(unsigned short dssd,unsigned short stripNo){
        fstripsmultx2[dssd]=fstripsmultx2[dssd]|((ULong64_t)0x1<<(stripNo-64));
    }
    void SetStripMultY1FlagMask(unsigned short dssd,unsigned short stripNo){
        fstripsmulty1[dssd]=fstripsmulty1[dssd]|((ULong64_t)0x1<<(stripNo-128));
    }
    void SetStripMultY2FlagMask(unsigned short dssd,unsigned short stripNo){
        fstripsmulty2[dssd]=fstripsmulty2[dssd]|((ULong64_t)0x1<<(stripNo-192));
    }



    //! Add a cluster
    void AddCluster(AIDACluster* cluster){
      cluster->SetClustersAdded(fnclusters);
      fclusters.push_back(cluster);
      int Z=(int) cluster->GetHitPositionZ();
      fnclustersz[Z]++;
      if (fnclustersz[Z]==1) fclustermultz++;
      fnclusters++;
    }

    //! Add more clusters
    void AddClusters(vector<AIDACluster*> clusters){
      fnclusters += clusters.size();
      for(vector<AIDACluster*>::iterator cluster=clusters.begin(); cluster!=clusters.end(); cluster++){
        fclusters.push_back(*cluster);
      }
    }

    //! Set all clusters
    void SetClusters(vector<AIDACluster*> clusters){
      fnclusters = clusters.size();
      fclusters = clusters;
    }

    //! Set number of cluster
    void SetNClusters(unsigned short ncluster){fnclusters = ncluster;}

    //! Set number of X clusters
    void SetNXClustersZ(int dssd, unsigned short ncluster){fnyclustersz[dssd] = ncluster;}

    //! Set number of Y clusters
    void SetNYClustersZ(int dssd, unsigned short ncluster){fnyclustersz[dssd] = ncluster;}

    //! Set number of dssd with at least 1 cluster
    void SetClusterMultZ(unsigned short clustermultz){
        fclustermultz = clustermultz;
    }

    //! Set number of cluster in dssd
    void SetNClusterZ(unsigned short nclusterz[NumDSSD]){
        memcpy(fnclustersz,nclusterz,NumDSSD*sizeof(unsigned short));
    }

    //! Set max Z (for ion event)
    void SetMaxZ(unsigned short maxz){fmaxz = maxz;}

    //! Set delta max Z (added July 29,2017)
    void SetDeltaMaxZ(unsigned short dmaxz){fdmaxz = dmaxz;}

    //! Set min distance between clusters
    void SetMinDistanceX(unsigned short dssd,double mindx){fmindx[dssd]=mindx;}
    void SetMinDistanceY(unsigned short dssd,double mindy){fmindy[dssd]=mindy;}


    //! Set event type
    void SetType(short type){ftype = type;}
    //! Set event type
    void SetTypeAdded(short type){ftype += type;}

    //! Set threshold

    //void SetThreshold(Double_t dssd_thr[NumDSSD][NumStrXY]){
    //    for (Int_t i=0;i<NumDSSD;i++){
    //        for (Int_t j=0;j<NumStrXY;j++){
    //            fdssd_thr[i][j] = dssd_thr[i][j];
    //        }
    //    }
    //}

    //! Set calibrationtable
    //void SetCalib(Double_t dssd_cal[NumDSSD][NumStrXY][2]){
    //    for (Int_t i=0;i<NumDSSD;i++){
    //        for (Int_t j=0;j<NumStrXY;j++){
    //            fdssd_cal[i][j][0] = dssd_cal[i][j][0];
    //            fdssd_cal[i][j][1] = dssd_cal[i][j][1];
    //        }
    //    }
    //}

    //! Returns timestamp
    unsigned long long GetTimestamp(){return faidats;}

    //! Revturn prev ion timestamp
    unsigned long long GetPrevIonTimestamp(){return ftsprevion;}

    //! Returns time width(window)
    unsigned long long GetTimeWindow(){return faidatw;}

    //! Returns the multiplicity of the event
    unsigned short GetMult(){return fmult;}

    //! Return the number of dssd with hit
    unsigned short GetZHitMult(){return fhitmultz;}

    //! Return the hits in z
    unsigned short* GetNHitZ(){return fnhitsz;}
    //! Returns the X strip multiplicity of the event
    unsigned short* GetMultXs(){return fmultx;}
    //! Returns the Y strip multiplicity of the event
    unsigned short* GetMultYs(){return fmulty;}

    //! Returns the X strip in dssd multiplicity of the event
    unsigned short GetMultX(short dssd){return fmultx[dssd];}
    //! Returns the Y strip in dssd multiplicity of the event
    unsigned short GetMultY(short dssd){return fmulty[dssd];}

    //! Returns the X strip in dssd sum energy of the event
    double GetSumX(short dssd){return fsumx[dssd];}
    //! Returns the Y strip in dssd sum energy of the event
    double GetSumY(short dssd){return fsumy[dssd];}

    //! Returns the number of dssd with cluster
    unsigned short GetClustersMultZ(){return fclustermultz;}

    unsigned short* GetNClustersZ(){return fnclustersz;}

    unsigned short GetNClustersZi(int dssd){return fnclustersz[dssd];}

    unsigned short GetNXClustersZi(int dssd){return fnxclustersz[dssd];}

    unsigned short GetNYClustersZi(int dssd){return fnyclustersz[dssd];}


    //! Returns the number of clusters
    unsigned short GetNClusters(){return fnclusters;}

    //! Returns the whole vector of clusters
    vector<AIDACluster*> GetClusters(){return fclusters;}
    //! Returns the hit number n
    AIDACluster* GetCluster(unsigned short n){return fclusters.at(n);}

    //! Clone Cluster
    AIDACluster* CloneCluster(unsigned short n){
        return (AIDACluster*) (fclusters.at(n)->Clone());
    }

    //! Return max Z (for ion event)
    unsigned short GetMaxZ(){return fmaxz;}

    //! Return delta max Z (july 27 added)
    unsigned short GetDeltaMaxZ(){return fdmaxz;}


    //! Returns the whole vector of hits
    vector<AIDAHit*> GetHits(){return fhits;}
    //! Returns the hit number n
    AIDAHit* GetHit(unsigned short n){return fhits.at(n);}

    //! Get threshold
    //Double_t GetThreshold(Int_t dssdno, Int_t stripno){return fdssd_thr[dssdno][stripno];}

    //! Get calibration table
    //Double_t GetCalibSlope(Int_t dssdno, Int_t stripno){return fdssd_cal[dssdno][stripno][0];}
    //! Get calibration table
    //Double_t GetCalibOffset(Int_t dssdno, Int_t stripno){return fdssd_cal[dssdno][stripno][1];}

    //! Get event type
    short GetType(){return ftype;}


    //! Get per silicon channel multi hit flag
    int GetStripMultFlag(unsigned short dssd){
        return fstripsmult[dssd];
    }

    unsigned long long GetStripMultX1Flag(unsigned short dssd){
        return fstripsmultx1[dssd];
    }
    unsigned long long GetStripMultX2Flag(unsigned short dssd){
        return fstripsmultx2[dssd];
    }
    unsigned long long GetStripMultY1Flag(unsigned short dssd){
        return fstripsmulty1[dssd];
    }
    unsigned long long GetStripMultY2Flag(unsigned short dssd){
        return fstripsmulty2[dssd];
    }

    //! get min distance between clusters
    double GetMinDistanceX(unsigned short dssd){
        return fmindx[dssd];
    }
    double GetMinDistanceY(unsigned short dssd){
        return fmindy[dssd];
    }

    //! Printing information
    void Print(Option_t *option = "") const {
      cout <<"timestamp " << faidats << endl;
        cout << "multiplicity " << fmult << endl;
      cout << "number of clusters " << fnclusters << endl;
      for(unsigned short i=0;i<NumDSSD;i++)
      cout << "DSSD " << i << " X strip multiplicity " << fmultx[i] << " Y strip multiplicity" << fmulty[i] << endl;
      for(unsigned short i=0;i<fhits.size();i++)
        fhits.at(i)->Print();
      for(unsigned short i=0;i<fclusters.size();i++)
        fclusters.at(i)->Print();
    }

    //! Get beta hit positions and fill in the cluster vector (clustering algorithms)
    //! Return true if there is at least 1 cluster identifed! otherwise return false

    //! Beta position new (with EX/EY condition)
    bool BetaGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);
    //! Get all Beta position new
    bool BetaGetPosAllNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);

    //! dev
    bool BetaGetPosAllNew2(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);

    //! Get all Beta position with position assignment using max. energy
    bool BetaGetPosAllNewMax(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);
    //! Ion position with clustering algorithm
    bool IonGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);
    bool IonGetPosAllNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);

    bool IonGetPosAllNew2(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);

    //! Get all ion position with position assignment using max. energy
    bool IonGetPosAllNewMax(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);


    //! only Y strips
    bool BetaGetPosYonly();
    bool ConstructEdiffRankingFlags();

  protected:
    //! aida time stamp (ealiest timestamp within event)
    unsigned long long faidats;

    //! the time window of the event (added July26)
    unsigned long long faidatw;

    //! timestamp of previous ion
    unsigned long long ftsprevion;

    //! per silicon channel multi hit flag
    int fstripsmult[NumDSSD];
    unsigned long long fstripsmultx1[NumDSSD];
    unsigned long long fstripsmultx2[NumDSSD];
    unsigned long long fstripsmulty1[NumDSSD];
    unsigned long long fstripsmulty2[NumDSSD];


    //! type of event: 0T beta  1 ion
    short ftype;
    //! total multiplicity
    unsigned short fmult;

    //! total multiplicity in z
    unsigned short fnhitsz[NumDSSD];

    //! number of dssd with at least 1 hit
    unsigned short fhitmultz;


    //! x multiplicity
    unsigned short fmultx[NumDSSD];
    //! y multiplicity
    unsigned short fmulty[NumDSSD];

    //! x sum all energy
    double fsumx[NumDSSD];
    //! y sum all energy
    double fsumy[NumDSSD];


    //! total clusters
    unsigned short fnclusters;

    //! total clusters
    unsigned short fnclustersz[NumDSSD];

    //! total X clusters
    unsigned short fnxclustersz[NumDSSD];
    //! total Y clusters
    unsigned short fnyclustersz[NumDSSD];

    //! max hit position
    unsigned short fmaxz;

    //! max hit position correction (delta_{maxz}) (added July 29,2017)
    unsigned short fdmaxz;

    //! number of dssd with at least 1 cluster
    unsigned short fclustermultz;

    //! Min distance between cluster
    double fmindx[NumDSSD];
    double fmindy[NumDSSD];

    //!Threshold table
    //Double_t fdssd_thr[NumDSSD][NumStrXY];

    //!Calibration table
    //Double_t fdssd_cal[NumDSSD][NumStrXY][2];

    //! vector with the hits
    vector<AIDAHit*> fhits;

    //! vector with the clusters
    vector<AIDACluster*> fclusters;

    /// \cond CLASSIMP
    ClassDef(AIDA,1);
    /// \endcond
};


#endif // AIDA_H
