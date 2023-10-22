#ifndef DATASTRUCTNEW_H
#define DATASTRUCTNEW_H
#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "Beam.h"
#include "TClonesArray.h"

class IonBeta : public TObject
{
public:
    IonBeta(){
        Clear();
    }
    virtual ~IonBeta(){}
    virtual void Clear(){
        for (size_t idx=0;idx<neuf.size();idx++){
            delete neuf[idx];
        }
        neuf.clear();
        for (size_t idx=0;idx<neub.size();idx++){
            delete neub[idx];
        }
        neub.clear();

        for (size_t idx=0;idx<clover.size();idx++){
            delete clover[idx];
        }
        clover.clear();
        for (size_t idx=0;idx<anc.size();idx++){
            delete anc[idx];
        }
        anc.clear();
        nneuf=0;
        nneub=0;
        nclover=0;
        nanc=0;
        nbeam=0;
        fid=0.;
        fevt=0;
        fts=0;
        fx=-9999.;
        fy=-9999.;
        fz=-1;
        fex=-9999.;
        fey=-9999.;
        fnx=0;
        fny=0;
        fncx=0;
        fncy=0;
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
        fminx = -11;
        fmaxx = -11;
        fminy = -11;
        fmaxy = -11;
        ftdiffmin= -99999;
        ftdiffmax= -99999;
        fdtionall= -9999.;
    }
    //! copy cluster
    virtual void Copy(IonBeta& obj){
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
        obj.SetXMinPos(fminx);
        obj.SetXMaxPos(fmaxx);
        obj.SetYMinPos(fminy);
        obj.SetYMaxPos(fmaxy);
        obj.SetTimeDifferenceMax(ftdiffmax);
        obj.SetTimeDifferenceMin(ftdiffmin);
    }
    void CopyFromAIDA(AIDASimpleStruct* obj){
        fid=obj->GetID();
        fevt=obj->GetEventNumber();
        fts=obj->GetTimestamp();
        fmult=obj->GetMultiplicity();
        fx=obj->GetHitPositionX();
        fy=obj->GetHitPositionY();
        fz=obj->GetHitPositionZ();
        fex=obj->GetXEnergy();
        fey=obj->GetYEnergy();
        fncx=obj->GetXClusterMultiplicity();
        fncy=obj->GetYClusterMultiplicity();
        fnx=obj->GetXMultiplicity();
        fny=obj->GetYMultiplicity();
        fnz=obj->GetZMultiplicity();
        fniz=obj->GetStripMultFlag();
        frflag=obj->GetRankingFlag();
        fdrflag=obj->GetEDiffRankingFlag();
        ftw=obj->GetTimeWidth();
        fdtion=obj->GetDtIon();
        fdz=obj->GetDZ();
        fmindx=obj->GetMinimumDistanceX();
        fmindy=obj->GetMinimumDistanceY();
        fsumexyrank=obj->GetSumEXYRank();
        fmaxelastdssd=obj->GetMaxELastDSSD();
        fxts=obj->GetXMinTimestamp();
        fyts=obj->GetYMinTimestamp();
        fminx=obj->GetMinHitPositionX();
        fmaxx=obj->GetMaxHitPositionX();
        fminy=obj->GetMinHitPositionY();
        fmaxy=obj->GetMaxHitPositionY();
        ftdiffmax=obj->GetTimeDifferenceMax();
        ftdiffmin=obj->GetTimeDifferenceMin();
    }

    void AddNeutronForward(BELENHit* neuhitin){
        neuf.push_back(neuhitin);
        nneuf++;
    }
    void AddNeutronBackward(BELENHit* neuhitin){
        neub.push_back(neuhitin);
        nneub++;
    }
    void AddAnc(BELENHit* anchitin){
        anc.push_back(anchitin);
        nanc++;
    }

    void AddClover(CloverHit* cloverhitin){
        clover.push_back(cloverhitin);
        nclover++;
    }

    BELENHit* GetNeutronForwardHit(unsigned short n){return neuf.at(n);}
    BELENHit* GetNeutronBackwardHit(unsigned short n){return neub.at(n);}
    CloverHit* GetCloverHit(unsigned short n){return clover.at(n);}
    BELENHit* GetAncHit(unsigned short n){return anc.at(n);}

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

    void SetEventNumber(unsigned int evt){fevt = evt;}

    //! Set the X strips sum energy
    void SetXEnergy(double xenergy){fex = xenergy;}
    //! Set the Y strips sum energy
    void SetYEnergy(double yenergy){fey = yenergy;}

    //! Set the X strips ealiest timestamp
    void SetXTimestamp(unsigned long long xts){fxts = xts;}
    //! Set the Y strips ealiest timestamp
    void SetYTimestamp(unsigned long long yts){fyts = yts;}

    //! Set time diffrernce between X-Y
    void SetTimeDifferenceMin(int tdiffmin){ftdiffmin = tdiffmin;}
    void SetTimeDifferenceMax(int tdiffmax){ftdiffmax = tdiffmax;}

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

    //! Set XY energy ratio ranking flag
   void SetRankingFlag(unsigned char rankingflag){frflag = rankingflag;}

   //! Set XY energy diffrence ranking flag
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

    //! Get the latest time since last ion-punching + trhough event and area
    void SetDtIonAll(double dtionall){fdtionall=dtionall;}

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


    //! Get the X strips ealiest timesttamp
    unsigned long long GetXTimestamp(){return fxts;}
    //! Get the Y strips ealiest timesttamp
    unsigned long long GetYTimestamp(){return fyts;}

    //! Get the time diffrerence of Y and X strips (min)
    int GetTimeDifferenceMin(){return ftdiffmin;}
    //! Get the time diffrerence of Y and X strips (max)
    int GetTimeDifferenceMax(){return ftdiffmax;}

    //! Get forward Neutron multiplicity
    unsigned short GetNeutronForwardMultipliticy(){return nneuf;}

    //! Get backward Neutron multiplicity
    unsigned short GetNeutronBackwardMultipliticy(){return nneub;}

    //! Get clover multiplicity
    unsigned short GetCloverMultipliticy(){return nclover;}

    //! Get beam multiplicity
    unsigned short GetBeamMultipliticy(){return nbeam;}

    //! Get Ancinary multiplicity
    unsigned short GetAncMultipliticy(){return nanc;}

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

    //! Get Max Energy of last layer
    double GetMaxELastDSSD(){return fmaxelastdssd;}

    //! Get the latest time since last ion-punching + trhough event and area
    double GetDtIonAll(){return fdtionall;}

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
    vector<BELENHit*> neuf;//forward correlation
    vector<BELENHit*> neub;//backward correlation
    vector<CloverHit*> clover;// correlation with anc detector
    vector<BELENHit*> anc;// correlation with anc detector
    unsigned short nneuf;
    unsigned short nneub;
    unsigned short nclover;
    unsigned short nanc;
    unsigned short nbeam;

    unsigned long long fts;

    //! id :ion or beta
    unsigned char fid;

    unsigned int fevt;

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

    //! the energy lab system
    unsigned long long fxts;
    unsigned long long fyts;

    //! time diffrence between y and x cluster
    int ftdiffmin;
    int ftdiffmax;


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

    //! energy diffrence ranking flag
    unsigned char fdrflag;

    unsigned short fsumexyrank;

    //! minimum cluster distance x
    double fmindx;
    //! minimum cluster distance y
    double fmindy;

    //! event time width in us
    double ftw;
    //! time distance with prev ion in us
    double fdtion;

    //! time distance with prev ion and puchching trhough events within a given area
    double fdtionall;

    //! z correction (if applicable)
    unsigned short fdz;

    //! max energy of the last layer (for light particle veto method)
    double fmaxelastdssd;

    /// \cond CLASSIMP
    ClassDef(IonBeta,1);
    /// \endcond
};

class IonBetaMult : public TObject
{
public:
    IonBetaMult(){
        //fnanc = 0;
        beam=new TreeData;
        Clear();
    }
    virtual ~IonBetaMult(){}
    virtual void Clear(){
        for (size_t idx=0;idx<neuf.size();idx++){
            delete neuf[idx];
        }
        neuf.clear();
        for (size_t idx=0;idx<neub.size();idx++){
            delete neub[idx];
        }
        neub.clear();

        beam->ts=0;
        beam->sts=0;
        beam->tof=-9999;
        beam->zet=-9999;
        beam->aoq=-9999;
        beam->f5x=-9999;
        beam->f11x=-9999;
        beam->f11y=-9999;
        beam->f11dt=-9999;
        beam->beta=-9999;

        for (size_t idx=0;idx<clover.size();idx++){
            delete clover[idx];
        }
        clover.clear();
        for (size_t idx=0;idx<anc.size();idx++){
            delete anc[idx];
        }
        anc.clear();
        nneuf=0;
        nneub=0;
        nclover=0;
        nanc=0;
        nbeam=0;

        fid=0.;
        fevt=0;
        fts=0;
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

        fxts=0;
        fyts=0;
        fminx = -11;
        fmaxx = -11;
        fminy = -11;
        fmaxy = -11;
        ftdiffmin= -99999;
        ftdiffmax= -99999;
        fdeltax = -99999.;
        fdeltay = -99999.;
    }

    virtual void Copy(IonBetaMult& obj){
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
        obj.SetXTimestamp(fxts);
        obj.SetYTimestamp(fyts);
        obj.SetXMinPos(fminx);
        obj.SetXMaxPos(fmaxx);
        obj.SetYMinPos(fminy);
        obj.SetYMaxPos(fmaxy);
        obj.SetTimeDifferenceMax(ftdiffmax);
        obj.SetTimeDifferenceMin(ftdiffmin);
        obj.SetDeltaX(fdeltax);
        obj.SetDeltaY(fdeltay);

        TreeData beamhitin;
        beamhitin.ts=beam->ts;
        beamhitin.ts=beam->sts;
        beamhitin.tof=beam->tof;
        beamhitin.zet=beam->zet;
        beamhitin.aoq=beam->aoq;
        beamhitin.f5x=beam->f5x;
        beamhitin.f11x=beam->f11x;
        beamhitin.f11y=beam->f11y;
        beamhitin.f11dt=beam->f11dt;
        beamhitin.beta=beam->beta;
        obj.AddBeam(beamhitin);

        for (unsigned short i=0;i<nneuf;i++){
            BELENHit* hit=neub.at(i);
            BELENHit* hitc=new BELENHit;
            hit->Copy(*hitc);
            obj.AddNeutronForward(hitc);
        }


        for (unsigned short i=0;i<nneub;i++){
            BELENHit* hit=neub.at(i);
            BELENHit* hitc=new BELENHit;
            hit->Copy(*hitc);
            obj.AddNeutronBackward(hitc);
        }

        for (unsigned short i=0;i<nclover;i++){
            CloverHit* hit=clover.at(i);
            CloverHit* hitc=new CloverHit;
            hit->Copy(*hitc);
            obj.AddClover(hitc);
        }

        for (unsigned short i=0;i<nanc;i++){
            BELENHit* hit=GetAncHit(i);
            BELENHit* hitc=new BELENHit;
            hit->Copy(*hitc);
            obj.AddAnc(hitc);
        }

    }

    //! copy cluster
    virtual void CopyWithBigRIPSOnly(IonBetaMult& obj){
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
        obj.SetXTimestamp(fxts);
        obj.SetYTimestamp(fyts);
        obj.SetXMinPos(fminx);
        obj.SetXMaxPos(fmaxx);
        obj.SetYMinPos(fminy);
        obj.SetYMaxPos(fmaxy);
        obj.SetTimeDifferenceMax(ftdiffmax);
        obj.SetTimeDifferenceMin(ftdiffmin);

        TreeData beamhitin;
        beamhitin.ts=beam->ts;
        beamhitin.ts=beam->sts;
        beamhitin.tof=beam->tof;
        beamhitin.zet=beam->zet;
        beamhitin.aoq=beam->aoq;
        beamhitin.f5x=beam->f5x;
        beamhitin.f11x=beam->f11x;
        beamhitin.f11y=beam->f11y;
        beamhitin.f11dt=beam->f11dt;
        beamhitin.beta=beam->beta;
        obj.AddBeam(beamhitin);        
    }
    void CopyFromAIDA(AIDASimpleStruct* obj){
        fid=obj->GetID();
        fevt=obj->GetEventNumber();
        fts=obj->GetTimestamp();
        fmult=obj->GetMultiplicity();
        fx=obj->GetHitPositionX();
        fy=obj->GetHitPositionY();
        fz=obj->GetHitPositionZ();
        fex=obj->GetXEnergy();
        fey=obj->GetYEnergy();
        fnx=obj->GetXMultiplicity();
        fny=obj->GetYMultiplicity();
        fncx=obj->GetXClusterMultiplicity();
        fncy=obj->GetYClusterMultiplicity();
        fnz=obj->GetZMultiplicity();
        fniz=obj->GetStripMultFlag();
        frflag=obj->GetRankingFlag();
        fdrflag=obj->GetEDiffRankingFlag();
        ftw=obj->GetTimeWidth();
        fdtion=obj->GetDtIon();
        fdz=obj->GetDZ();
        fmindx=obj->GetMinimumDistanceX();
        fmindy=obj->GetMinimumDistanceY();
        fsumexyrank=obj->GetSumEXYRank();
        fxts=obj->GetXMinTimestamp();
        fyts=obj->GetYMinTimestamp();
        fminx=obj->GetMinHitPositionX();
        fmaxx=obj->GetMaxHitPositionX();
        fminy=obj->GetMinHitPositionY();
        fmaxy=obj->GetMaxHitPositionY();
        ftdiffmax=obj->GetTimeDifferenceMax();
        ftdiffmin=obj->GetTimeDifferenceMin();
    }

    void AddNeutronForward(BELENHit* neuhitin){
        neuf.push_back(neuhitin);
        nneuf++;
    }
    void AddNeutronBackward(BELENHit* neuhitin){
        neub.push_back(neuhitin);
        nneub++;
    }
    void AddAnc(BELENHit* anchitin){
        anc.push_back(anchitin);
        nanc++;
    }
    void AddBeam(TreeData& beamhitin){
        beam->ts=beamhitin.ts;
        beam->sts=beamhitin.ts;
        beam->tof=beamhitin.tof;
        beam->zet=beamhitin.zet;
        beam->aoq=beamhitin.aoq;
        beam->f5x=beamhitin.f5x;
        beam->f11x=beamhitin.f11x;
        beam->f11y=beamhitin.f11y;
        beam->f11dt=beamhitin.f11dt;
        beam->beta=beamhitin.beta;
        nbeam++;
    }
    void AddClover(CloverHit* cloverhitin){
        clover.push_back(cloverhitin);
        nclover++;
    }

    BELENHit* GetNeutronForwardHit(unsigned short n){return neuf.at(n);}
    BELENHit* GetNeutronBackwardHit(unsigned short n){return neub.at(n);}
    CloverHit* GetCloverHit(unsigned short n){return clover.at(n);}
    BELENHit* GetAncHit(unsigned short n){return anc.at(n);}
    TreeData* GetBeamHit(){return beam;}


    unsigned short GetNNeutronForwardHit(){return nneuf;}
    unsigned short GetNNeutronBackwardHit(){return nneub;}
    unsigned short GetNCloverHit(){return nclover;}
    unsigned short GetNAncHit(){return nanc;}


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

    void SetEventNumber(unsigned int evt){fevt = evt;}

    //! Set the X strips sum energy
    void SetXEnergy(double xenergy){fex = xenergy;}
    //! Set the Y strips sum energy
    void SetYEnergy(double yenergy){fey = yenergy;}

    //! Set the X strips ealiest timestamp
    void SetXTimestamp(unsigned long long xts){fxts = xts;}
    //! Set the Y strips ealiest timestamp
    void SetYTimestamp(unsigned long long yts){fyts = yts;}

    //! Set time diffrernce between X-Y
    void SetTimeDifferenceMin(int tdiffmin){ftdiffmin = tdiffmin;}
    void SetTimeDifferenceMax(int tdiffmax){ftdiffmax = tdiffmax;}

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

    //! Set XY energy ratio ranking flag
    void SetRankingFlag(unsigned char rankingflag){frflag = rankingflag;}


    //! Set XY energy diffrence ranking flag
    void SetEDiffRankingFlag(unsigned char rankingflag){fdrflag = rankingflag;}

    //! Set EX+EY rank
    void SetSumEXYRank(unsigned short sumexyrank){fsumexyrank = sumexyrank;}

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

    //! Set distance between clusters
    void SetDeltaX(double deltax){fdeltax=deltax;}
    void SetDeltaY(double deltay){fdeltay=deltay;}


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

    //! Get the X strips ealiest timesttamp
    unsigned long long GetXTimestamp(){return fxts;}
    //! Get the Y strips ealiest timesttamp
    unsigned long long GetYTimestamp(){return fyts;}

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

    //! Get Sum energy ranking
    unsigned short GetSumEXYRank(){return fsumexyrank;}

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

    //! Get distance between clusters
    double GetDeltaX(){return fdeltax;}
    double GetDeltaY(){return fdeltay;}


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
    vector<BELENHit*> neuf;//forward correlation
    vector<BELENHit*> neub;//backward correlation
    vector<CloverHit*> clover;// correlation with anc detector
    vector<BELENHit*> anc;// correlation with anc detector
    TreeData* beam;//correlation with beam in bigrip
    unsigned short nneuf;
    unsigned short nneub;
    unsigned short nclover;
    unsigned short nanc;
    unsigned short nbeam;


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

    //! ealiest timestamp in 1 cluster
    unsigned long long fxts;
    unsigned long long fyts;

    //! time diffrence between y and x cluster
    int ftdiffmin;
    int ftdiffmax;

    //! xy hit multiplicity
    unsigned short fnx;
    unsigned short fny;\
    //! xy cluster multiplicity
    unsigned short fncx;
    unsigned short fncy;

    //! number of clusters
    unsigned short fnz;
    //! number of event with mult>1 for 1 dssd
    unsigned short fniz;
    //! ranking flag
    unsigned char frflag;

    //! ediff ranking flag
    unsigned char fdrflag;

    unsigned short fsumexyrank;

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



    //! global description of delta xy
    //! negative value: beta and ion area is overlap the value is the overlaping range
    //! positive value: beta and ion area is not overlap, the value is the minimum distance
    //! defined as ratio between cluster center and sum of the width of implant-decay cluster
    double fdeltax;
    double fdeltay;


    /// \cond CLASSIMP
    ClassDef(IonBetaMult,1);
    /// \endcond
};

class AncNeutron : public TObject
{
public:
    AncNeutron(){
        Clear();
    }
    virtual ~AncNeutron(){}
    virtual void Clear(){
        for (size_t idx=0;idx<neu.size();idx++){
            delete neu[idx];
        }
        neu.clear();
        anc->Clear();
        nneu=0;
    }

    void AddNeutron(BELENHit* neuhitin){
        neu.push_back(neuhitin);
        nneu++;
    }

    //! Get forward Neutron multiplicity
    unsigned short GetNeutronMultipliticy(){return nneu;}
    BELENHit* GetAnc(){return anc;}
    BELENHit* GetNeutronHit(unsigned short n){return neu.at(n);}

    //! Printing information
    void Print(Option_t *option = "") const {

    }

protected:
    BELENHit* anc;
    vector<BELENHit*> neu;
    unsigned short nneu;

    /// \cond CLASSIMP
    ClassDef(AncNeutron,1);
    /// \endcond
};



#endif //  DATASTRUCTNEW_H
