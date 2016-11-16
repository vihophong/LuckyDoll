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

      memset(fsumx,0,sizeof(fsumx));
      memset(fsumy,0,sizeof(fsumy));


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
      if (hit->GetXY() < NumStrX) {
          fmultx[z]++;
          fsumx[z]+=hit->GetEnergy();
      }else {
          fmulty[z]++;
          fsumy[z]+=hit->GetEnergy();
      }
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

    //! Returns the X strip in dssd sum energy of the event
    double GetSumX(short dssd){return fsumx[dssd];}
    //! Returns the Y strip in dssd sum energy of the event
    double GetSumY(short dssd){return fsumy[dssd];}


    unsigned short* GetNClustersZ(){return fnclustersz;}

    unsigned short GetNClustersZi(int dssd){return fnclustersz[dssd];}

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

bool AIDA::BetaGetPos(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    int maxmult=1;
    Double_t dssdH_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdH_E_Y[NumDSSD][NumStrY][maxmult];
    Long64_t dssdH_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdH_T_Y[NumDSSD][NumStrY][maxmult];

    Long64_t dssdF_T_X[NumDSSD][NumStrX][maxmult];
    Long64_t dssdF_T_Y[NumDSSD][NumStrY][maxmult];

    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdH_E_X[i][j][k]=0;
                dssdH_E_Y[i][j][k]=0;
                dssdH_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
                dssdF_T_X[i][j][k]=0;
                dssdF_T_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        unsigned long long time = hit->GetTimestamp();
        unsigned long long fastime = hit->GetFastTimestamp();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) {
                dssdH_E_X[z][xy][mult_strip_x[z][xy]] = energy;
                dssdH_T_X[z][xy][mult_strip_x[z][xy]] = time;
                dssdF_T_X[z][xy][mult_strip_x[z][xy]] = fastime;
            }
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) {
                dssdH_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
                dssdH_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = time;
                dssdF_T_Y[z][xy-128][mult_strip_y[z][xy-128]] = fastime;
            }
            mult_strip_y[z][xy-128]++;
        }
    }

    //! Old code

    Int_t maxStrInCluster = 3;
    Int_t maxNposBeta = 3;
    Int_t maxNCluster = 4;

    for(Int_t z=0; z<NumDSSD; z++){
        Int_t mult_z = 0;
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > > beta_dssd_pre;
        //vector<pair<E-corr,pair<pair<PosX,pair<EX,NClusterX> >,pair<PosY,pair<EY,NClusterY> > >  > >
        vector<pair<Double_t,pair<pair<Double_t,pair<Double_t,Int_t> >,pair<Double_t,pair<Double_t,Int_t> > >  > >::iterator ibeta_dssd_pre;

        Double_t posX=-1;
        Double_t posY=-1;

        Int_t nClusterX=0;
        Int_t nClusterY=0;
        Int_t nStripsX=0;
        Double_t E_X=0;
        Double_t E_X_ch=0;

        for(Int_t x=0; x<NumStrX; x++){
            //cout << "x-strip found" << endl;
            if (dssdH_E_X[z][x][0]>0) {
                E_X+=dssdH_E_X[z][x][0];
                E_X_ch+=x*dssdH_E_X[z][x][0];
                nStripsX++;
            }
            if ((dssdH_E_X[z][x][0]<=0||x==NumStrX-1)&&nStripsX>0) {//end of X cluster
                if (nStripsX<=maxStrInCluster){ //if cluster of less than 3 strips
                    posX=E_X_ch/E_X;
                    //another loop for Y strips
                    nClusterY=0;
                    Int_t nStripsY=0;
                    Double_t E_Y=0;
                    Double_t E_Y_ch=0;
                    for(Int_t y=0; y<NumStrY; y++){
                        //cout << "x-strip found" << endl;
                        if (dssdH_E_Y[z][y][0]>0) {
                            E_Y+=(Double_t) dssdH_E_Y[z][y][0];
                            E_Y_ch+=(Double_t) y*dssdH_E_Y[z][y][0];
                            nStripsY++;
                        }
                        if ((dssdH_E_Y[z][y][0]<=0||y==NumStrY-1)&&nStripsY>0) {//end of X cluster
                            if (nStripsY<=maxStrInCluster){ //if cluster of less than 3 strips
                                posY=E_Y_ch/E_Y;
                                beta_dssd_pre.push_back(make_pair((E_Y/E_X-1)*(E_Y/E_X-1),make_pair(make_pair(posX,make_pair(E_X,nStripsX)),make_pair(posY,make_pair(E_Y,nStripsY)))));
                            }
                            E_Y=0;
                            E_Y_ch=0;
                            nStripsY=0;
                            nClusterY++;
                        }
                    }
                    //ok!
                }
                E_X=0;
                E_X_ch=0;
                nStripsX=0;
                nClusterX++;
            }
        }

        //Additional step to filter all possible combination:
        if (nClusterX>0&&nClusterX<maxNCluster&&nClusterY>0&&nClusterY<maxNCluster){
            //Record number of cluster here
            vector <Double_t> xindex;
            vector <Double_t> yindex;
            sort(beta_dssd_pre.begin(),beta_dssd_pre.end());
            for (ibeta_dssd_pre=beta_dssd_pre.begin();ibeta_dssd_pre<beta_dssd_pre.end();++ibeta_dssd_pre){
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end()))||
                        ((corr_cut>0)&&(ibeta_dssd_pre->first<corr_cut*corr_cut)&&(find(xindex.begin(),xindex.end(),ibeta_dssd_pre->second.first.first)==xindex.end())&&(find(yindex.begin(),yindex.end(),ibeta_dssd_pre->second.second.first)==yindex.end())))
                {
                    //beta_dssd.push_back(make_pair(make_pair(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.first.second),make_pair(ibeta_dssd_pre->second.second.first,ibeta_dssd_pre->second.second.second)));
                    if (mult_z<maxNposBeta && ibeta_dssd_pre->second.first.second.first>sumexcut[z] && ibeta_dssd_pre->second.second.second.first>sumeycut[z]){
                        AIDACluster *cluster=new AIDACluster;
                        cluster->SetHitPosition(ibeta_dssd_pre->second.first.first,ibeta_dssd_pre->second.second.first,z);
                        cluster->SetXEnergy(ibeta_dssd_pre->second.first.second.first);
                        cluster->SetYEnergy(ibeta_dssd_pre->second.second.second.first);
                        cluster->SetXMult(ibeta_dssd_pre->second.first.second.second);
                        cluster->SetYMult(ibeta_dssd_pre->second.second.second.second);
                        //if (ibeta_dssd_pre->second.first.second.second!=ibeta_dssd_pre->second.second.second.second) cout<<"eurica!"<<cluster->GetXMultiplicity()<<"-"<<cluster->GetYMultiplicity()<<endl;

                        //Timing
                        int hitx = (int) round(ibeta_dssd_pre->second.first.first);
                        int hity = (int) round(ibeta_dssd_pre->second.second.first);
                        if (dssdH_T_X[z][hitx][0]>0&&dssdH_T_Y[z][hity][0]>0){

                            //! if there is available fast time stamp
                            if (dssdH_T_Y[z][hitx][0]==0||dssdH_T_Y[z][hity][0]==0){
                                cluster->SetFastTimestamp(dssdH_T_Y[z][hitx][0] + dssdH_T_Y[z][hity][0]);
                            }
                            //! Take the ealiest time stamp
                            if (dssdH_T_X[z][hitx][0]>dssdH_T_Y[z][hity][0]) {
                                cluster->SetTimestamp(dssdH_T_Y[z][hity][0]);
                                cluster->SetFastTimestamp(dssdF_T_Y[z][hity][0]);
                            }else {
                                cluster->SetTimestamp(dssdH_T_X[z][hitx][0]);
                                cluster->SetFastTimestamp(dssdF_T_X[z][hitx][0]);
                            }

                        }else{
                            cout<<__PRETTY_FUNCTION__<<" something wrong!"<<endl;
                        }
                        this->AddCluster(cluster);
                    }
                    xindex.push_back(ibeta_dssd_pre->second.first.first);
                    yindex.push_back(ibeta_dssd_pre->second.second.first);
                    mult_z++;
                }
            }
        }//end of additional step...

    }//end of loop on all dssd
    if(this->GetNClusters()>0) return true;
    return false;
}

bool AIDA::IonGetPos()
{
    int maxmult=1;
    Double_t dssdL_E_X[NumDSSD][NumStrX][maxmult];
    Double_t dssdL_E_Y[NumDSSD][NumStrY][maxmult];
    int mult_strip_x[NumDSSD][NumStrX];
    int mult_strip_y[NumDSSD][NumStrX];
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrX;j++){
            mult_strip_x[i][j]=0;
            mult_strip_y[i][j]=0;
            for (int k=0;k<maxmult;k++){
                dssdL_E_X[i][j][k]=0;
                dssdL_E_Y[i][j][k]=0;
            }
        }
    }
    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        int xy = hit->GetXY();
        double energy = hit->GetEnergy();
        if (xy<128){
            if (mult_strip_x[z][xy]<maxmult) dssdL_E_X[z][xy][mult_strip_x[z][xy]] = energy;
            mult_strip_x[z][xy]++;
        }else{
            if (mult_strip_y[z][xy-128]<maxmult) dssdL_E_Y[z][xy-128][mult_strip_y[z][xy-128]] = energy;
            mult_strip_y[z][xy-128]++;
        }
    }
    //!Old code


    //get pos by weighting method
    Double_t Ion_posX;
    Double_t Ion_posY;
    Bool_t is_X_hit=false;
    Bool_t is_Y_hit=false;
    Int_t my_ion_mult_x,my_ion_mult_y;
    Int_t maxZ=-1;

    for(Int_t z=0; z<NumDSSD; z++){
        is_X_hit=false;
        is_Y_hit=false;
        my_ion_mult_x=0;
        my_ion_mult_y=0;
        Double_t sum_ex=0;
        Double_t sum_weightx=0;
        Double_t sum_ey=0;
        Double_t sum_weighty=0;
        for(Int_t x=0; x<NumStrX; x++){
            if(dssdL_E_X[z][x][0]>0){
                is_X_hit=true;
                sum_ex+=dssdL_E_X[z][x][0];
                sum_weightx+=x*dssdL_E_X[z][x][0];
                my_ion_mult_x++;
            }
        }
        if (is_X_hit){
            for(Int_t y=0; y<NumStrY; y++){
                if(dssdL_E_Y[z][y][0]>0){
                   is_Y_hit=true;
                   sum_ey+=dssdL_E_Y[z][y][0];
                   sum_weighty+=y*dssdL_E_Y[z][y][0];
                   my_ion_mult_y++;
               }
            }
            if(is_Y_hit){
                Ion_posX=sum_weightx/sum_ex;
                Ion_posY=sum_weighty/sum_ey;
                //! only one cluster
                AIDACluster *cluster=new AIDACluster;
                cluster->SetXEnergy(sum_ex);
                cluster->SetYEnergy(sum_ey);
                cluster->SetHitPosition(Ion_posX,Ion_posY,z);
                cluster->SetXMult(my_ion_mult_x);
                cluster->SetYMult(my_ion_mult_y);
                //! take ealiest time stamp
                cluster->SetTimestamp(fhits.at(0)->GetTimestamp());
                this->AddCluster(cluster);
                maxZ=z;
            }
        }
    }
    this->SetMaxZ(maxZ);
    if (maxZ!=-1) return true;
    return false;
}

