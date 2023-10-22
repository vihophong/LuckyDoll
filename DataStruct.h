#ifndef DATASTRUCT_H
#define DATASTRUCT_H
#include "AIDA.h"
#include "BELEN.h"
#include "Clover.h"
#include "Beam.h"

class Ancillary : public TObject
{
public:
    Ancillary(){
        fnf11r = 0;
        fnf11l = 0;
        fnvetobot = 0;
        fnvetotop = 0;
        fnvetodown = 0;
        fe_f11l = 0;
        fe_f11r = 0;
        fe_vetobot = 0;
        fe_vetotop = 0;
        fe_vetodown = 0;
        ft_f11l = 0;
        ft_f11r = 0;
        ft_vetobot = 0;
        ft_vetotop = 0;
        ft_vetodown = 0;
        fnneutron = 0;

        fndebot = 0;
        fndetop = 0;
        fe_detop = 0;
        fe_debot = 0;
        ft_detop = 0;
        ft_debot = 0;
    }
    virtual ~Ancillary(){}
    virtual void Clear(){
        fnf11r = 0;
        fnf11l = 0;
        fnvetobot = 0;
        fnvetotop = 0;
        fnvetodown = 0;
        fe_f11l = 0;
        fe_f11r = 0;
        fe_vetobot = 0;
        fe_vetotop = 0;
        fe_vetodown = 0;
        ft_f11l = 0;
        ft_f11r = 0;
        ft_vetobot = 0;
        ft_vetotop = 0;
        ft_vetodown = 0;
        fnneutron = 0;
        fngamma = 0;

        fndebot = 0;
        fndetop = 0;
        fe_detop = 0;
        fe_debot = 0;
        ft_detop = 0;
        ft_debot = 0;

        //! Dealocating memory
        for (size_t idx=0;idx<fNeutrons.size();idx++){
            delete fNeutrons[idx];
        }
        fNeutrons.clear();
        for (size_t idx=0;idx<fClovers.size();idx++){
            delete fClovers[idx];
        }
        fClovers.clear();
    }
    //! Copy object
    virtual void Copy(Ancillary& obj){
        obj.SetNF11L(fnf11l);
        obj.SetNF11R(fnf11r);
        obj.SetNVetoTop(fnvetotop);
        obj.SetNVetoBot(fnvetobot);
        obj.SetNVetoDown(fnvetodown);

        obj.SetEF11L(fe_f11l);
        obj.SetEF11R(fe_f11r);
        obj.SetEVetoTop(fe_vetotop);
        obj.SetEVetoBot(fe_vetobot);
        obj.SetEVetoDown(fe_vetodown);

        obj.SetTF11L(ft_f11l);
        obj.SetTF11R(ft_f11r);
        obj.SetTVetoTop(ft_vetotop);
        obj.SetTVetoBot(ft_vetobot);
        obj.SetTVetoDown(ft_vetodown);

        obj.SetNdETop(fndetop);
        obj.SetNdEBot(fndebot);
        obj.SetTdETop(ft_detop);
        obj.SetTdEBot(ft_debot);
        obj.SetEdETop(fe_detop);
        obj.SetEdEBot(fe_debot);
        for (size_t idx=0;idx<fNeutrons.size();idx++){
            BELENHit* hitneutron=new BELENHit;
            fNeutrons[idx]->Copy(*hitneutron);
            obj.AddNeutronHit(hitneutron);
        }
        for (size_t idx=0;idx<fClovers.size();idx++){
            CloverHit* hitclover=new CloverHit;
            fClovers[idx]->Copy(*hitclover);
            obj.AddCloverHit(hitclover);
        }


    }

    void SetNF11L(unsigned char nf11l){fnf11l = nf11l;}
    void SetNF11R(unsigned char nf11r){fnf11r = nf11r;}
    void SetNVetoTop(unsigned char nvetotop){fnvetotop = nvetotop;}
    void SetNVetoBot(unsigned char nvetobot){fnvetobot = nvetobot;}
    void SetNVetoDown(unsigned char nvetodown){fnvetodown = nvetodown;}

    void SetNdETop(unsigned char ndetop){fndetop = ndetop;}
    void SetNdEBot(unsigned char ndebot){fndebot = ndebot;}

    unsigned char GetNF11L(){return fnf11l;}
    unsigned char GetNF11R(){return fnf11r;}
    unsigned char GetNVetoTop(){return fnvetotop;}
    unsigned char GetNVetoBot(){return fnvetobot;}
    unsigned char GetNVetoDown(){return fnvetodown;}



    double GetEF11L(){return fe_f11l;}
    double GetEF11R(){return fe_f11r;}
    double GetEVetoTop(){return fe_vetotop;}
    double GetEVetoBot(){return fe_vetobot;}
    double GetEVetoDown(){return fe_vetodown;}

    unsigned long long GetTF11L(){return ft_f11l;}
    unsigned long long GetTF11R(){return ft_f11r;}
    unsigned long long GetTVetoTop(){return ft_vetotop;}
    unsigned long long GetTVetoBot(){return ft_vetobot;}
    unsigned long long GetTVetoDown(){return ft_vetodown;}

    double GetEdETop(){return fe_detop;}
    double GetTdETop(){return ft_detop;}
    double GetEdEBot(){return fe_debot;}
    double GetTdEBot(){return ft_debot;}

    void SetEF11L(double ef11l){fe_f11l = ef11l;}
    void SetEF11R(double ef11r){fe_f11r = ef11r;}
    void SetEVetoTop(double evetotop){fe_vetotop = evetotop;}
    void SetEVetoBot(double evetobot){fe_vetobot = evetobot;}
    void SetEVetoDown(double evetodown){fe_vetodown = evetodown;}

    void SetTF11L(unsigned long long tf11l){ft_f11l = tf11l;}
    void SetTF11R(unsigned long long tf11r){ft_f11r = tf11r;}
    void SetTVetoTop(unsigned long long tvetotop){ft_vetotop = tvetotop;}
    void SetTVetoBot(unsigned long long tvetobot){ft_vetobot = tvetobot;}
    void SetTVetoDown(unsigned long long tvetodown){ft_vetodown = tvetodown;}

    void SetEdETop(double edetop){fe_detop = edetop;}
    void SetTdETop(double tdetop){ft_detop = tdetop;}
    void SetEdEBot(double edebot){fe_debot = edebot;}
    void SetTdEBot(double tdebot){ft_debot = tdebot;}

    unsigned long long GetTimeStamp(){return ft_f11l;}

    //! Add Neutron hits
    void AddNeutronHit(BELENHit* hit){
        fNeutrons.push_back(hit);
        fnneutron++;
    }
    //! Add Clover hits
    void AddCloverHit(CloverHit* hit){
        fClovers.push_back(hit);
        fngamma++;
    }

private:
    unsigned char fnf11l;
    unsigned char fnf11r;
    unsigned char fnvetotop;
    unsigned char fnvetobot;
    unsigned char fnvetodown;

    unsigned char fndetop;
    unsigned char fndebot;

    double fe_f11r;
    double fe_f11l;
    double fe_vetotop;
    double fe_vetobot;
    double fe_vetodown;

    double fe_detop;
    double fe_debot;


    unsigned long long ft_f11r;
    unsigned long long ft_f11l;
    unsigned long long ft_vetotop;
    unsigned long long ft_vetobot;
    unsigned long long ft_vetodown;

    unsigned long long ft_detop;
    unsigned long long ft_debot;

    //! BELEN hits
    unsigned short fnneutron;
    vector<BELENHit*> fNeutrons;

    //! Clover hits
    unsigned short fngamma;
    vector<CloverHit*> fClovers;

    /// \cond CLASSIMP
    ClassDef(Ancillary,1);
    /// \endcond
};


class Implant : public TObject
{
public:
    Implant(){
        ftimestamp = 0;
        fnbeam = 0;
        fnneutron = 0;
        fngamma = 0;
        //fnanc = 0;
        fBeam = new Beam;
        //fIon = new AIDACluster;
        fIon = new AIDA;
        fF11Beam = new Ancillary;
    }
    virtual ~Implant(){}
    virtual void Clear(){
        ftimestamp = 0;
        fnbeam = 0;
        fnneutron = 0;
        fngamma = 0;
        //fnanc = 0;
        fIon->Clear();
        fBeam->Clear();
        //! Dealocating memory
        for (size_t idx=0;idx<fNeutrons.size();idx++){
            delete fNeutrons[idx];
        }
        fNeutrons.clear();
        for (size_t idx=0;idx<fClovers.size();idx++){
            delete fClovers[idx];
        }
        fClovers.clear();
        //for (size_t idx=0;idx<fAncs.size();idx++){
        //    delete fAncs[idx];
        //}
        //fAncs.clear();
        fF11Beam->Clear();
    }
    //! Set time stamp of this beta event
    void SetTimeStamp(unsigned long long ts){ftimestamp = ts;}
    //! Get time stamp of this beta event
    unsigned long long GetTimeStamp(){return ftimestamp;}
    //! Get Cluster
    //AIDACluster* GetIon(){return fIon;}
    AIDA* GetIon(){return fIon;}
    //! Get Beam
    Beam* GetBeam(){return fBeam;}
    //! Get F11 Beam
    Ancillary* GetF11Beam(){return fF11Beam;}
    //! Set number of correlated beam (0 or 1)
    void SetNBeam(unsigned short nbeam){fnbeam=nbeam;}

    //! Add Neutron hits
    void AddNeutronHit(BELENHit* hit){
        fNeutrons.push_back(hit);
        fnneutron++;
    }
    //! Add Clover hits
    void AddGammaHit(CloverHit* hit){
        fClovers.push_back(hit);
        fngamma++;
    }
    //! Add Anc hits
    //void AddAncHit(BELENHit* hit){
    //    fAncs.push_back(hit);
    //    fnanc++;
    //}

protected:
    //! timestamp
    unsigned long long ftimestamp;
    unsigned short fnbeam;
    unsigned short fnneutron;
    unsigned short fngamma;
    //unsigned short fnanc;

    //! AIDA implantation
    AIDA* fIon;

    //! Bigrips Beam
    Beam* fBeam;
    //! BELEN hits
    vector<BELENHit*> fNeutrons;
    //! Clover hits
    vector<CloverHit*> fClovers;
    //! Anc hits
    //vector<BELENHit*> fAncs;
    //! Beam hits (coincidence anc hits)
    Ancillary* fF11Beam;

    /// \cond CLASSIMP
    ClassDef(Implant,1);
    /// \endcond
};

class Beta : public TObject
{
public:
    Beta(){
        ftimestamp = 0;
        fngamma = 0;
        //fnanc = 0;
        fmult = 0;
        fnclusters = 0;
        fmaxz = -1;
        fclustermultz = 0;
        fhitmultz = 0;
        memset(fmultx,0,sizeof(fmultx));
        memset(fmulty,0,sizeof(fmulty));
        memset(fnhitsz,0,sizeof(fnhitsz));
        memset(fnclustersz,0,sizeof(fnclustersz));
        fBetas = new AIDACluster;
        fF11Beam = new Ancillary;
    }
    virtual ~Beta(){}
    virtual void Clear(){
        ftimestamp = 0;
        fngamma = 0;
        //fnanc = 0;

        fmult = 0;
        fnclusters = 0;
        fmaxz = -1;
        fclustermultz = 0;
        fhitmultz = 0;
        memset(fmultx,0,sizeof(fmultx));
        memset(fmulty,0,sizeof(fmulty));
        memset(fnhitsz,0,sizeof(fnhitsz));
        memset(fnclustersz,0,sizeof(fnclustersz));

        fBetas->Clear();
        //! Dealocating memory
        for (size_t idx=0;idx<fClovers.size();idx++){
            delete fClovers[idx];
        }
        fClovers.clear();
        //for (size_t idx=0;idx<fAncs.size();idx++){
        //    delete fAncs[idx];
        //}
        //fAncs.clear();
        fF11Beam->Clear();
    }
    //! Set time stamp of this beta event
    void SetTimeStamp(unsigned long long ts){ftimestamp = ts;}
    //! Get time stamp of this beta event
    unsigned long long GetTimeStamp(){return ftimestamp;}

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


    //! Get beta
    AIDACluster* GetBeta(){return fBetas;}
    //! Get F11 Beam
    Ancillary* GetF11Beam(){return fF11Beam;}
    //! Get Clover Beam
    CloverHit* GetCloverHit(unsigned int i){return fClovers.at(i);}
    unsigned short GetCloverMult(){return fngamma;}
    //! Get Clover mult

    //! Add Clover hits
    void AddGammaHit(CloverHit* hit){
        fClovers.push_back(hit);
        fngamma++;
    }

    //! Add Anc hits
    //void AddAncHit(BELENHit* hit){
    //    fAncs.push_back(hit);
    //    fnanc++;
    //}


    //!common stuff for light particle rejection later
    //! all multiplicity
    unsigned short fmult;
    //! number of dssd with a hit
    unsigned short fhitmultz;
    //! total multiplicity in z
    unsigned short fnhitsz[NumDSSD];
    //! x multiplicity
    unsigned short fmultx[NumDSSD];
    //! y multiplicity
    unsigned short fmulty[NumDSSD];
    //! max z
    unsigned short fmaxz;
    //! total clusters
    unsigned short fnclusters;
    //! total clusters in z
    unsigned short fnclustersz[NumDSSD];
    //! number of dssd with a cluster
    unsigned short fclustermultz;


protected:
    unsigned long long ftimestamp;
    unsigned short fngamma;
    //unsigned short fnanc;

    /*
    //!common stuff for light particle rejection later
    //! all multiplicity
    unsigned short fmult;
    //! number of dssd with a hit
    unsigned short fhitmultz;
    //! total multiplicity in z
    unsigned short fnhitsz[NumDSSD];
    //! x multiplicity
    unsigned short fmultx[NumDSSD];
    //! y multiplicity
    unsigned short fmulty[NumDSSD];
    //! max z
    unsigned short fmaxz;
    //! total clusters
    unsigned short fnclusters;
    //! total clusters in z
    unsigned short fnclustersz[NumDSSD];
    //! number of dssd with a cluster
    unsigned short fclustermultz;
    */


    //! AIDA cluster
    AIDACluster* fBetas;
    //! Clover hits
    vector<CloverHit*> fClovers;
    //! Anc hits
    //vector<BELENHit*> fAncs;

    //! Beam hits (coincidence anc hits)
    Ancillary* fF11Beam;

    /// \cond CLASSIMP
    ClassDef(Beta,1);
    /// \endcond
};

class Neutron : public TObject
{
public:
    Neutron(){
        ftimestamp = 0;
        fngamma = 0;
        //fnanc = 0;
        fNeutrons=new BELENHit;
        fF11Beam = new Ancillary;
    }
    virtual ~Neutron(){}
    virtual void Clear(){
        ftimestamp = 0;
        fngamma = 0;
        //fnanc = 0;
        fNeutrons->Clear();
        //! Dealocating memory
        for (size_t idx=0;idx<fClovers.size();idx++){
            delete fClovers[idx];
        }
        fClovers.clear();
        //for (size_t idx=0;idx<fAncs.size();idx++){
        //    delete fAncs[idx];
        //}
        //fAncs.clear();
        fF11Beam->Clear();

    }
    //! Copy object
    virtual void Copy(Neutron& obj){
        obj.SetTimeStamp(ftimestamp);
        obj.SetNGamma(fngamma);
        //obj.SetNAnc(fnanc);
        fNeutrons->Copy(*obj.GetNeutron());
        for(vector<CloverHit*>::iterator clover=fClovers.begin(); clover!=fClovers.end(); clover++){
            CloverHit* cloneclover=new CloverHit;
            CloverHit* originclover = *clover;
            originclover->Copy(*cloneclover);
            obj.AddGammaHit(cloneclover);
        }
        //for(vector<BELENHit*>::iterator anc=fAncs.begin(); anc!=fAncs.end(); anc++){
        //    BELENHit* cloneanc=new BELENHit;
        //    BELENHit* originanc = *anc;
        //    originanc->Copy(*cloneanc);
        //    obj.AddAncHit(cloneanc);
        //}
        fF11Beam->Copy(*obj.GetF11Beam());
    }
    //! Set time stamp of this beta event
    void SetTimeStamp(unsigned long long ts){ftimestamp = ts;}
    //! Get time stamp of this beta event
    unsigned long long GetTimeStamp(){return ftimestamp;}

    void SetNGamma(unsigned short ngamma){fngamma = ngamma;}
    //void SetNAnc(unsigned short nanc){fnanc = nanc;}

    //! Get neutron
    BELENHit* GetNeutron(){return fNeutrons;}

    //! Get F11 Beam
    Ancillary* GetF11Beam(){return fF11Beam;}

    //! Add Clover hits
    void AddGammaHit(CloverHit* hit){
        fClovers.push_back(hit);
        fngamma++;
    }
    //! Add Anc hits
    //void AddAncHit(BELENHit* hit){
    //    fAncs.push_back(hit);
    //    fnanc++;
    //}
protected:
    unsigned long long ftimestamp;
    unsigned short fngamma;
    //unsigned short fnanc;
    //! AIDA cluster
    BELENHit* fNeutrons;
    //! Clover hits
    vector<CloverHit*> fClovers;
    //! Ancillary hits
    //vector<BELENHit*> fAncs;
    //! Beam hits (coincidence anc hits)
    Ancillary* fF11Beam;
    /// \cond CLASSIMP
    ClassDef(Neutron,1);
    /// \endcond
};

class Neutrons : public TObject
{
public:
    Neutrons(){
        fnneutron = 0;
    }
    virtual ~Neutrons(){}
    virtual void Clear(){
        fnneutron = 0;
        //! Dealocating memory
        for (size_t idx=0;idx<fNeutrons.size();idx++){
            delete fNeutrons[idx];
        }
        fNeutrons.clear();
    }
    //! Add Anc hits
    void AddNeutron(Neutron* hit){
        fNeutrons.push_back(hit);
        fnneutron++;
    }
    unsigned int GetMult(){return fnneutron;}
    Neutron* GetNeutron(int i){return fNeutrons.at(i);}
private:
    unsigned int fnneutron;
    vector<Neutron*> fNeutrons;
    /// \cond CLASSIMP
    ClassDef(Neutrons,1);
    /// \endcond
};





#endif //  DATASTRUCT_H
