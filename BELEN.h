#ifndef BELEN_H
#define BELEN_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "BELENdefs.h"

using namespace std;


class BELENHit : public TObject{
public:
    BELENHit(){
        Clear();
    }
    BELENHit(Double_t posx, Double_t posy, Double_t posz, short id, unsigned long long ts, int adc, int en, unsigned short hitsadded)
    {
        fid = id;
        fts = ts;
        fadc = adc;
        fen = en;
        fhe3pos.SetXYZ(posx,posy,posz);
        fhitsadded = hitsadded;
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
    void SetPos(Double_t x, Double_t y, Double_t z){fhe3pos.SetXYZ(x,y,z);}
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
    TVector3 GetPosition(){return fhe3pos;}

    //! Get current hits
    unsigned short GetHitsAdded(){return fhitsadded;}

    void Clear(){
        fid = -1;
        fhe3pos.SetXYZ(-1,-1,-1);
        fts = 0;
        fadc = -1;
        fen = -1;
        fhitsadded = 0;
    }

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "ID " << fid;
      cout << "\tX pos " << fhe3pos.X();
      cout << "\tY pos " << fhe3pos.Y();
      cout << "\tZ pos " << fhe3pos.Z();
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
    TVector3 fhe3pos;
    //! ID number of 3He counter
    short fid;
    //! ADC value
    int fadc;
    //! Energy calibrated value
    double fen;
    //! timestamp value
    unsigned long long fts;

    /// \cond CLASSIMP
    ClassDef(BELENHit,1);
    /// \endcond
    ///
};

class BELEN : public TObject
{
public:
    //! default constructor
    BELEN(){
        Clear();
    }
    //! clear BELEN info
    void Clear(){
        fmult = 0;
        //! Dealocating memory
        for (size_t idx=0;idx<fhits.size();idx++){
            delete fhits[idx];
        }
        fhits.clear();
    }

    //! Set time stamp
    void SetTimestamp(unsigned long long ts){fbelents = ts;}
    //! Set Multiplicity;
    void SetMult(unsigned short mult) {fmult = mult;}

    //! Add more hits
    void AddHits(vector<BELENHit*> hits){
      fmult += hits.size();
      for(vector<BELENHit*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
        //set hit add here!
        //
        fhits.push_back(*hit);
      }
    }

    //! Add a hit
    void AddHit(BELENHit* hit){
      //!newly added
      hit->SetHitsAdded(fmult);
      fhits.push_back(hit);
      fmult++;
    }

    //! Set all hits
    void SetHits(vector<BELENHit*> hits){
      fmult = hits.size();
      fhits = hits;
    }

    //! Returns timestamp
    unsigned long long GetTimestamp(){return fbelents;}
    //! Returns the multiplicity of the event
    unsigned short GetMult(){return fmult;}


    //! Returns the whole vector of hits
    vector<BELENHit*> GetHits(){return fhits;}
    //! Returns the hit number n
    BELENHit* GetHit(unsigned short n){return fhits.at(n);}


    void Print(Option_t *option = "") const {
        cout <<"timestamp " << fbelents << endl;
        cout << "multiplicity " << fmult << endl;
    }
protected:
    //! total multiplicity
    unsigned short fmult;
    //! the ealiest time stamp found
    unsigned long long  fbelents;
    //! vector with the hits
    vector<BELENHit*> fhits;

    /// \cond CLASSIMP
    ClassDef(BELEN,1);
    /// \endcond
};

#endif // BELLEN_H
