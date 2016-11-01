#ifndef CLOVER_H
#define CLOVER_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "Cloverdefs.h"

using namespace std;

class CloverHit : public TObject
{
public:
    CloverHit(){
        Clear();
    }
    void Clear(){
        fid = -1;
        fpos.SetXYZ(-1,-1,-1);
        fts = 0;
        fadc = -1;
        fen = -1;
        fhitsadded = 0;
    }
    //! Set the energy
    void SetEnergy(double energy){fen = energy;}

    //! Set the raw ADC value
    void SetADC(int adc){fadc = adc;}

    //! Set the counter ID
    void SetID(short id){fid = id;}

    //! Set the timestamp
    void SetTimestamp(unsigned long long int ts){fts = ts;}

    //! Set the He3 position
    void SetPos(Double_t x, Double_t y, Double_t z){fpos.SetXYZ(x,y,z);}
    //! Set current hits
    void SetHitsAdded(unsigned short hitsadded){fhitsadded = hitsadded;}


    //! Get the ID
    short GetID(){return fid;}
    //! Get the energy
    double GetEnergy(){return fen;}
    //! Get the timestamp
    unsigned long long int GetTimestamp(){return fts;}
    //! Get the raw ADC value
    int GetADC(){return fadc;}

    //! Get 3He position
    TVector3 GetPosition(){return fpos;}

    //! Get current hits
    unsigned short GetHitsAdded(){return fhitsadded;}

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "ID " << fid;
      cout << "\tX pos " << fpos.X();
      cout << "\tY pos " << fpos.Y();
      cout << "\tZ pos " << fpos.Z();
      cout << "\tadc " << fadc;
      cout << "\tenergy " << fen;
      cout << "\ttimestamp " << fts;
      cout << "\thits added " << fhitsadded << endl;
      return;
    }

protected:
    //! current hits
    unsigned short fhitsadded;
    //! Position of 3He counter
    TVector3 fpos;
    //! ID number of 3He counter
    short fid;
    //! ADC value
    int fadc;
    //! Energy calibrated value
    double fen;
    //! timestamp value
    unsigned long long fts;

    /// \cond CLASSIMP
    ClassDef(CloverHit,1);
    /// \endcond
    ///
};

class Clover : public TObject
{
public:
    Clover(){
        Clear();
    }
    //! Clear the music information
    void Clear(Option_t *option = ""){
      fmult = 0;
      fmultAB = 0;
      //! Dealocating memory
      for (size_t idx=0;idx<fhits.size();idx++){
          delete fhits[idx];
      }
      fhits.clear();
      for (size_t idx=0;idx<fhitsAB.size();idx++){
          delete fhitsAB[idx];
      }
      fhitsAB.clear();
    }
    //! Add a hit
    void AddHit(CloverHit* hit){
      hit->SetHitsAdded(fmult);
      fhits.push_back(hit);
      fmult++;
    }
    //! Add a hit after addback
    void AddHitAB(CloverHit* hit){
      fhitsAB.push_back(hit);
      fmultAB++;
    }
    //! Set all hits
    void SetHits(vector<CloverHit*> hits){
      fmult = hits.size();
      fhits = hits;
    }
    //! Set all hits after addback
    void SetABHits(vector<CloverHit*> hits){
      fmultAB = hits.size();
      fhitsAB = hits;
    }
    //! Add more hits
    void AddHits(vector<CloverHit*> hits){
      fmult += hits.size();
      for(vector<CloverHit*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
          //set hit add here!
          //
          fhits.push_back(*hit);
      }
    }


    //! Returns the multiplicity of the event
    unsigned short GetMult(){return fmult;}
    //! Returns the whole vector of hits
    vector<CloverHit*> GetHits(){return fhits;}
    //! Returns the hit number n
    CloverHit* GetHit(unsigned short n){return fhits.at(n);}
    //! Returns the multiplicity of the event after addback
    int GetMultAB(){return fmultAB;}
    //! Returns the whole vector of hits after addback
    vector<CloverHit*> GetHitsAB(){return fhitsAB;}
    //! Returns the hit number n after addback
    CloverHit* GetHitAB(int n){return fhitsAB.at(n);}

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "multiplicity " << fmult << " event" << endl;
      for(unsigned short i=0;i<fhits.size();i++)
        fhits.at(i)->Print();
      cout << "after addback multiplicity " << fmultAB << endl;
      for(unsigned short i=0;i<fhitsAB.size();i++)
        fhitsAB.at(i)->Print();
    }
protected:
    //! multiplicity
    unsigned short fmult;
    //! vector with the hits
    vector<CloverHit*> fhits;
    //! multiplicity after addback
    unsigned short fmultAB;
    //! vector with the hits after addback
    vector<CloverHit*> fhitsAB;


    /// \cond CLASSIMP
    ClassDef(Clover,1);
    /// \endcond
    ///
};

#endif // CLOVER_H
