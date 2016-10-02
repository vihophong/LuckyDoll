#ifndef AIDA_H
#define AIDA_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "AIDAdefs.h"
using namespace std;


class AIDAHit : public TObject {
public:
  //! default constructor
  AIDAHit(){
    Clear();
  }
  virtual ~AIDAHit(){}
  //! constructor with individual values
  AIDAHit(short range, short id, int xy, int z, double en, int adc,unsigned short hitsadded, unsigned long long int ts){
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
  void Clear(Option_t *option = ""){
      fid = 0;
      fxy = 0;
      fz = 0;
      fadc = -9999;
      fen = -9999.;
      fhitsadded = 0;
      fts = 0;
      ffastts = 0;
      frange = 0;
  }

  //! Set energy range (low gain 0 or high gain 1)
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
  //! the energy lab system
  double fen;
  //! the raw adc value
  int fadc;
  //! the timestamp
  unsigned long long fts;
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
              double xenergy,double yenergy,unsigned short clustersadded, unsigned long long int ts)
  {
    fpos.SetXYZ(x,y,z);
    fsumenx=xenergy;
    fsumeny=yenergy;
    fnx=multx;
    fny=multy;
    fclustersadded = clustersadded;
    fcts = ts;
  }
  //! Clear the music information
  void Clear(Option_t *option = ""){
      fpos.SetXYZ(-1,-1,-1);
      fsumenx=-9999.;
      fsumeny=-9999.;
      fcts=0;
      fclustersadded=0;
  }
  //! Set XY
  void SetHitPosition(double x,double y,double z){fpos.SetXYZ(x,y,z);}

  //! Set the X strips sum energy
  void SetXEnergy(double xenergy){fsumenx = xenergy;}
  //! Set the Y strips sum energy
  void SetYEnergy(double yenergy){fsumeny = yenergy;}
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


  TVector3 GetHitPosition(){return fpos;}
  Double_t GetHitPositionX(){return fpos.X();}
  Double_t GetHitPositionY(){return fpos.Y();}
  Double_t GetHitPositionZ(){return fpos.Z();}

  //! Get the X strips sum energy
  double GetXEnergy(){return fsumenx;}
  //! Get the Y strips sum energy
  double GetYEnergy(){return fsumeny;}

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
    return;
  }

protected:
  //! translate into the DSSDs coordinator
  TVector3 fpos;

  //! store the hit number
  vector<short*> fhitsno;

  //! the energy lab system
  double fsumenx;
  double fsumeny;

  //! number of hit in one cluster
  unsigned short fnx;
  unsigned short fny;

  //! the ealiest timestamp
  unsigned long long fcts;

  //! the ealiest fast timestamp
  unsigned long long fcfastts;

  //! current hits
  unsigned short fclustersadded;

  /// \cond CLASSIMP
  ClassDef(AIDACluster,1);
  /// \endcond
};

class AIDA : public TObject
{
public:
    //! default constructor
    AIDA(){
      for (Int_t i=0;i<NumDSSD;i++){
          for (Int_t j=0;j<NumStrXY;j++){
              fdssd_thr[i][j] = 0.;
              fdssd_cal[i][j][0] = 0.;
              fdssd_cal[i][j][1] = 1.;
          }
      }
      Clear();
    }
    virtual ~AIDA(){}
    //! Clear the AIDA information
    void Clear(Option_t *option = ""){

      fmult = 0;
      fnclusters = 0;
      fmaxz = -1;
      /*
      for (int i=0;i<NumDSSD;i++){
          fmultx[i]=0;
          fmulty[i]=0;
      }
      */
      memset(fmultx,0,sizeof(fmultx));
      memset(fmulty,0,sizeof(fmulty));

      memset(fnclustersz,0,sizeof(fnclustersz));


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

    void SetTimestamp(unsigned long long ts){faidats = ts;}

    //! Add a hit
    void AddHit(AIDAHit* hit){
      //!newly added
      hit->SetHitsAdded(fmult);
      fhits.push_back(hit);
      //! newly added
      Int_t z=(Int_t)hit->GetZ();
      if (hit->GetXY() < NumStrX) fmultx[z]++;
      else fmulty[z]++;
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

    //! Set the X strip multiplicity of the event
    void SetMultX(unsigned short multx[NumDSSD]){
        memcpy(fmultx,multx,NumDSSD*sizeof(unsigned short));
    }

    //! Set the Y strip multiplicity of the event
    void SetMultY(unsigned short multy[NumDSSD]){
        memcpy(fmulty,multy,NumDSSD*sizeof(unsigned short));
    }


    //! Add a cluster
    void AddCluster(AIDACluster* cluster){
      cluster->SetClustersAdded(fnclusters);
      fclusters.push_back(cluster);
      int Z=(int) cluster->GetHitPositionZ();
      fnclustersz[Z]++;
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
    //! Set max Z (for ion event)
    void SetMaxZ(unsigned short maxz){fmaxz = maxz;}
    //! Set event type
    void SetType(short type){ftype = type;}
    //! Set event type
    void SetTypeAdded(short type){ftype += type;}

    //! Set threshold
    void SetThreshold(Double_t dssd_thr[NumDSSD][NumStrXY]){
        for (Int_t i=0;i<NumDSSD;i++){
            for (Int_t j=0;j<NumStrXY;j++){
                fdssd_thr[i][j] = dssd_thr[i][j];
            }
        }
    }

    //! Set calibrationtable
    void SetCalib(Double_t dssd_cal[NumDSSD][NumStrXY][2]){
        for (Int_t i=0;i<NumDSSD;i++){
            for (Int_t j=0;j<NumStrXY;j++){
                fdssd_cal[i][j][0] = dssd_cal[i][j][0];
                fdssd_cal[i][j][1] = dssd_cal[i][j][1];
            }
        }
    }

    //! Returns timestamp
    unsigned long long GetTimestamp(){return faidats;}
    //! Returns the multiplicity of the event
    unsigned short GetMult(){return fmult;}

    //! Returns the X strip multiplicity of the event
    unsigned short* GetMultXs(){return fmultx;}
    //! Returns the Y strip multiplicity of the event
    unsigned short* GetMultYs(){return fmulty;}

    //! Returns the X strip in dssd multiplicity of the event
    unsigned short GetMultX(short dssd){return fmultx[dssd];}
    //! Returns the Y strip in dssd multiplicity of the event
    unsigned short GetMultY(short dssd){return fmulty[dssd];}


    unsigned short* GetNClustersZ(){return fnclustersz;}

    //! Returns the number of clusters
    unsigned short GetNClusters(){return fnclusters;}

    //! Returns the whole vector of clusters
    vector<AIDACluster*> GetClusters(){return fclusters;}
    //! Returns the hit number n
    AIDACluster* GetCluster(unsigned short n){return fclusters.at(n);}

    //! Return max Z (for ion event)
    unsigned short GetMaxZ(){return fmaxz;}

    //! Returns the whole vector of hits
    vector<AIDAHit*> GetHits(){return fhits;}
    //! Returns the hit number n
    AIDAHit* GetHit(unsigned short n){return fhits.at(n);}

    //! Get threshold
    Double_t GetThreshold(Int_t dssdno, Int_t stripno){return fdssd_thr[dssdno][stripno];}

    //! Get calibration table
    Double_t GetCalibSlope(Int_t dssdno, Int_t stripno){return fdssd_cal[dssdno][stripno][0];}
    //! Get calibration table
    Double_t GetCalibOffset(Int_t dssdno, Int_t stripno){return fdssd_cal[dssdno][stripno][1];}

    //! Get event type
    short GetType(){return ftype;}

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

    //! Beta position
    bool BetaGetPos(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);
    //! Ion position
    bool IonGetPos();

  protected:
    unsigned long long faidats;
    //! type of event: 0 beta  1 ion
    short ftype;
    //! total multiplicity
    unsigned short fmult;

    //! x multiplicity
    unsigned short fmultx[NumDSSD];
    //! y multiplicity
    unsigned short fmulty[NumDSSD];

    //! total clusters
    unsigned short fnclusters;

    //! total clusters
    unsigned short fnclustersz[NumDSSD];

    //! max hit position
    unsigned short fmaxz;

    //!Threshold table
    Double_t fdssd_thr[NumDSSD][NumStrXY];

    //!Calibration table
    Double_t fdssd_cal[NumDSSD][NumStrXY][2];


    //! vector with the hits
    vector<AIDAHit*> fhits;
    //! vector with the clusters
    vector<AIDACluster*> fclusters;

    /// \cond CLASSIMP
    ClassDef(AIDA,1);
    /// \endcond
};


#endif // AIDA_H
