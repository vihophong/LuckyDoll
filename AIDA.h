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
  virtual void Clear(){
      fpos.SetXYZ(-1,-1,-1);
      fsumenx=-9999.;
      fsumeny=-9999.;
      fcts=0;
      fclustersadded=0;
  }

  //! copy cluster
  virtual void Copy(AIDACluster& obj){
      obj.SetTimestamp(fcts);
      obj.SetHitPosition(fpos.X(),fpos.Y(),fpos.Z());
      obj.SetXEnergy(fsumenx);
      obj.SetYEnergy(fsumeny);
      obj.SetXMult(fnx);
      obj.SetYMult(fny);
      obj.SetFastTimestamp(fcfastts);
      obj.SetClustersAdded(fclustersadded);
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
  //vector<short*> fhitsno;

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
      fmult = 0;
      fnclusters = 0;
      fmaxz = -1;
      fclustermultz = 0;
      fhitmultz = 0;
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
          clonecluster->SetTimestamp(clonecluster->GetTimestamp()*ClockResolution);

	  //if (!(clonecluster->GetHitPositionZ()==1&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()<64))
	  //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>125))	    
	  //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>50&&clonecluster->GetHitPositionY()<75))
	  //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()<64&&clonecluster->GetHitPositionY()>98&&clonecluster->GetHitPositionY()<111))
	  //if (!(clonecluster->GetHitPositionZ()==0&&clonecluster->GetHitPositionX()>126))
          obj.AddCluster(clonecluster);
        }
        obj.SetTimestamp(faidats*ClockResolution);
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
    //! Beta position new
    bool BetaGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[]);
    //! Ion position
    bool IonGetPos();

  protected:
    //! aida time stamp (ealiest timestamp within event)
    unsigned long long faidats;
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

    //! max hit position
    unsigned short fmaxz;

    //! number of dssd with at least 1 cluster
    unsigned short fclustermultz;


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
