#include "Merger.h"

Merger::Merger()
{
    fIonPidTWup = 0;
    fIonPidTWlow = 50000;
    fIonNeutronTWup = 500000;
    fIonNeutronTWlow = 400000;
    fIonGammaTWup = 5000;
    fIonGammaTWlow = 70000;
    fIonAncTWup = 20000;
    fIonAncTWlow = 100000;

    fNeuGammaTWup = 50000;
    fNeuGammaTWlow = 500000;
    fNeuAncTWup = 50000;
    fNeuAncTWlow = 500000;

    fBetaGammaTWup = 40000;
    fBetaGammaTWlow = 50000;
    fBetaAncTWup = 40000;
    fBetaAncTWlow = 50000;

    fF11LRTWup = 400;
    fF11LRTWlow = 400;
    fF11LVetoTWup = 400;
    fF11LVetoTWlow  = 400;
    fF11LVetoDownup  = 1200;
    fF11LVetoDownlow  = -600;
}

Merger::~Merger()
{
    delete fh1;
}

void Merger::Init()
{
    ftreeImplant = 0;
    ftreeBeta = 0;
    ftreeNeutron = 0;
    flocalbeta = new Beta;
    flocalimp = new Implant;
    flocalneutron = new Neutron;
    flocalancillary = new Ancillary;

    faidaIon = NULL;
    faidaBeta = NULL;
    fbigrips = NULL;
    fclover = NULL;
    fneutron = NULL;
    fanc = NULL;

    fnentriesAIDAIon = 0;
    fnentriesAIDABeta = 0;
    fnentriesBigrips = 0;
    fnentriesGamma = 0;
    fnentriesNeutron = 0;
    fnentriesAnc = 0;

    //! init aida
    fAidaFile = new TFile(finputAida);
    fAidaFile->GetObject("beta",ftrAIDABeta);
    fAidaFile->GetObject("ion",ftrAIDAIon);
    ftrAIDABeta->SetBranchAddress("aida",&faidaBeta);
    ftrAIDAIon->SetBranchAddress("aida",&faidaIon);
    fnentriesAIDAIon = ftrAIDAIon->GetEntries();
    fnentriesAIDABeta = ftrAIDABeta->GetEntries();
    cout<<"Reading "<<fnentriesAIDABeta<<" betas and "
       <<fnentriesAIDAIon<<" ions in AIDA tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrAIDABeta->GetEvent(i);
        ftrAIDAIon->GetEvent(i);
        cout<<"betats "<<faidaBeta->GetTimestamp()<<"- ionts "<<faidaIon->GetTimestamp()<<endl;
    }

    //! init briken
    fBrikenFile = new TFile(finputBriken);
    fBrikenFile->GetObject("neutron",ftrNeutron);
    fBrikenFile->GetObject("gamma",ftrGamma);
    fBrikenFile->GetObject("anc",ftrAnc);
    ftrNeutron->SetBranchAddress("neutron",&fneutron);
    ftrGamma->SetBranchAddress("gamma",&fclover);
    ftrAnc->SetBranchAddress("anc",&fanc);
    fnentriesNeutron = ftrNeutron->GetEntries();
    fnentriesGamma = ftrGamma->GetEntries();
    fnentriesAnc = ftrAnc->GetEntries();
    cout<<"Reading "<<fnentriesNeutron<<" neutrons, "
       <<fnentriesGamma<<" gammas and "<<fnentriesAnc<<" anc hits in AIDA tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrNeutron->GetEvent(i);
        ftrGamma->GetEvent(i);
        ftrAnc->GetEvent(i);
        cout<<"neutronts "<<fneutron->GetTimestamp()<<"- gammats "<<fclover->GetTimestamp()<<"- ancts"<<fanc->GetTimestamp()<<endl;
    }

    //! init bigrips
    fBigripsFile = new TFile(finputBigrips);
    cout<<finputBigrips<<endl;
    fBigripsFile->GetObject("tree",ftrBigrips);
    ftrBigrips->SetBranchAddress("bigrips",&fbigrips);
    fnentriesBigrips = ftrBigrips->GetEntries();

    cout<<"Reading "<<fnentriesBigrips<<" bigrips items"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrBigrips->GetEvent(i);
        cout<<"bigripsts "<<fbigrips->GetTimestamp()<<endl;
    }
    cout<<"\n\t\t\t****************************************\n"<<endl;

    fh1=new TH1F("h1","h1",500,0,100e3);
}


void Merger::InitBRIKEN()
{
    ftreeImplant = 0;
    ftreeBeta = 0;
    ftreeNeutron = 0;
    flocalneutron = new Neutron;

    faidaIon = NULL;
    faidaBeta = NULL;
    fbigrips = NULL;
    fclover = NULL;
    fneutron = NULL;
    fanc = NULL;

    fnentriesAIDAIon = 0;
    fnentriesAIDABeta = 0;
    fnentriesBigrips = 0;
    fnentriesGamma = 0;
    fnentriesNeutron = 0;
    fnentriesAnc = 0;


    //! init briken
    fBrikenFile = new TFile(finputBriken);
    fBrikenFile->GetObject("neutron",ftrNeutron);
    fBrikenFile->GetObject("gamma",ftrGamma);
    fBrikenFile->GetObject("anc",ftrAnc);
    ftrNeutron->SetBranchAddress("neutron",&fneutron);
    ftrGamma->SetBranchAddress("gamma",&fclover);
    ftrAnc->SetBranchAddress("anc",&fanc);
    fnentriesNeutron = ftrNeutron->GetEntries();
    fnentriesGamma = ftrGamma->GetEntries();
    fnentriesAnc = ftrAnc->GetEntries();
    cout<<"Reading "<<fnentriesNeutron<<" neutrons, "
       <<fnentriesGamma<<" gammas and "<<fnentriesAnc<<" anc hits in AIDA tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrNeutron->GetEvent(i);
        ftrGamma->GetEvent(i);
        ftrAnc->GetEvent(i);
        cout<<"neutronts "<<fneutron->GetTimestamp()<<"- gammats "<<fclover->GetTimestamp()<<"- ancts"<<fanc->GetTimestamp()<<endl;
    }
    fh1=new TH1F("h1","h1",500,0,100e3);
}
void Merger::BookTreeBRIKEN(TTree* treeAnc){
    ftreeAncillary = treeAnc;
    ftreeAncillary->Branch("beam",&flocalancillary);
    ftreeAncillary->BranchRef();
}


/*
void Merger::InitTemp()
{
    flocalancillary = new Ancillary;
    fanc = NULL;
    fnentriesAnc = 0;
    //! init briken
    fBrikenFile = new TFile(finputBriken);
    fBrikenFile->GetObject("anc",ftrAnc);
    ftrAnc->SetBranchAddress("anc",&fanc);
    fnentriesAnc = ftrAnc->GetEntries();
    fh1=new TH1F("h1","h1",500,0,100e3);
}

void Merger::BookTreeTemp(TTree* treeAnc){
    ftreeAncillary = treeAnc;
    ftreeAncillary->Branch("beam",&flocalancillary);
}
*/

void Merger::BookTree(TTree* treeImplant, TTree* treeBeta, TTree* treeNeutron, TTree* treeAnc){
    ftreeImplant = treeImplant;
    ftreeBeta = treeBeta;
    ftreeNeutron = treeNeutron;
    ftreeAncillary = treeAnc;

    ftreeImplant->Branch("implant",&flocalimp);
    ftreeImplant->BranchRef();
    ftreeBeta->Branch("beta",&flocalbeta);
    ftreeBeta->BranchRef();
    ftreeNeutron->Branch("neutron",&flocalneutron);
    ftreeNeutron->BranchRef();
    ftreeAncillary->Branch("beam",&flocalancillary);
    ftreeAncillary->BranchRef();
}

void Merger::ReadAIDA()
{
    //! read ion
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesAIDAIon;jentry++){
        ftrAIDAIon->GetEvent(jentry);
        unsigned short lastclusterID = faidaIon->GetNClusters()-1;
        if (faidaIon->GetCluster(lastclusterID)->GetHitPositionZ()==faidaIon->GetMaxZ()){
            faidaIonMap.insert(make_pair(faidaIon->GetCluster(lastclusterID)->GetTimestamp(),jentry));
        }else{
            cout<<"something wrong!"<<endl;
        }
    }
    cout<<"Finished reading ion  ts table with "<<faidaIonMap.size()<<" rows"<<endl;
    //! read beta
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesAIDABeta;jentry++){
        ftrAIDABeta->GetEvent(jentry);
        for (unsigned short i = 0;i<faidaBeta->GetNClusters();i++){
            faidaBetaMap.insert(make_pair(faidaBeta->GetCluster(i)->GetTimestamp(),make_pair(jentry,i)));
        }
    }
    cout<<"Finished reading beta ts table with "<<faidaBetaMap.size()<<" rows"<<endl;
}

void Merger::ReadBigrips()
{
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesBigrips;jentry++){
        ftrBigrips->GetEvent(jentry);
        fbigripsMap.insert(make_pair(fbigrips->GetTimestamp(),jentry));
    }
    cout<<"Finished reading bigrips ts table with "<<fbigripsMap.size()<<" rows"<<endl;
}

void Merger::ReadBRIKEN()
{
    //! read neutron
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesNeutron;jentry++){
        ftrNeutron->GetEvent(jentry);
        fhe3Map.insert(make_pair(fneutron->GetTimestamp(),jentry));
    }
    cout<<"Finished reading neutron  ts table with "<<fhe3Map.size()<<" rows"<<endl;

    //! read gamma
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesGamma;jentry++){
        ftrGamma->GetEvent(jentry);
        fcloverMap.insert(make_pair(fclover->GetTimestamp(),jentry));
    }
    cout<<"Finished reading gamma  ts table with "<<fcloverMap.size()<<" rows"<<endl;

    //! read anc
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesAnc;jentry++){
        ftrAnc->GetEvent(jentry);
        //fancMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==1) fF11RMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==2) fF11LMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==1) fVetoTopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==2) fVetoBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==4) fVetoDownMap.insert(make_pair(fanc->GetTimestamp(),jentry));
    }
    //cout<<"Finished reading anc  ts table with "<<fancMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11R  ts table with "<<fF11RMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11L  ts table with "<<fF11LMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Top ts table with "<<fVetoTopMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Bottom ts table with "<<fVetoBotMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Down ts table with "<<fVetoDownMap.size()<<" rows"<<endl;
}

void Merger::ReadTemp()
{
    //! read anc
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesAnc;jentry++){
        ftrAnc->GetEvent(jentry);
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==1) fF11RMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==2) fF11LMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==1) fVetoTopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==2) fVetoBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==4) fVetoDownMap.insert(make_pair(fanc->GetTimestamp(),jentry));
    }
    cout<<"Finished reading F11R  ts table with "<<fF11RMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11L  ts table with "<<fF11LMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Top ts table with "<<fVetoTopMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Bottom ts table with "<<fVetoBotMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Down ts table with "<<fVetoDownMap.size()<<" rows"<<endl;
}

void Merger::DoMergeImp()
{
    unsigned int ncorrwBR = 0;
    unsigned int ncorrwNeutron = 0;
    unsigned int ncorrwGamma = 0;
    unsigned int ncorrwAnc = 0;
    unsigned int ncorrwF11R = 0;
    unsigned int ncorrwF11L = 0;
    unsigned int ncorrwVetoTop = 0;
    unsigned int ncorrwVetoBot = 0;
    unsigned int ncorrwVetoDown = 0;


    for (faidaIonMap_it=faidaIonMap.begin();faidaIonMap_it!=faidaIonMap.end();faidaIonMap_it++){
        unsigned long long ts  = faidaIonMap_it->first;
        unsigned int entry = faidaIonMap_it->second;
        ftrAIDAIon->GetEvent(entry);

        flocalimp->Clear();

        //! fill ion first
        unsigned short lastclusterID = faidaIon->GetNClusters()-1;
        if (faidaIon->GetCluster(lastclusterID)->GetHitPositionZ()==faidaIon->GetMaxZ()){
            faidaIon->GetCluster(lastclusterID)->Copy(*flocalimp->GetIon());
        }else{
            cout<<"something wrong!"<<endl;
        }
        flocalimp->SetTimeStamp(ts);
        //! correlated event with bigrips
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        fbigripsMap_it = fbigripsMap.lower_bound(ts1);
        while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
            corrts = (Long64_t) fbigripsMap_it->first;
            correntry = fbigripsMap_it->second;
            if (corrts!=check_time){
                ftrBigrips->GetEvent(correntry);
                fbigrips->Copy(*flocalimp->GetBeam());               
                check_time=corrts;
                flocalimp->SetNBeam(1);
                fh1->Fill(-(Long64_t)flocalimp->GetBeam()->GetTimestamp()+(Long64_t)flocalimp->GetIon()->GetTimestamp());
                 ncorrwBR++;
                break;
            }
            fbigripsMap_it++;
        }

        //! correlated event with neutron
        ts1 = (Long64_t)ts - (Long64_t)fIonNeutronTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonNeutronTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            correntry = fhe3Map_it->second;
            if (corrts!=check_time){
                ftrNeutron->GetEvent(correntry);
                check_time=corrts;
                BELENHit* hit=new BELENHit;
                fneutron->Copy(*hit);
                flocalimp->AddNeutronHit(hit);
                ncorr++;
            }
            fhe3Map_it++;
        }
        if (ncorr>0) ncorrwNeutron++;

        //! correlated event with gamma
        ts1 = (Long64_t)ts - (Long64_t)fIonGammaTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonGammaTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fcloverMap_it = fcloverMap.lower_bound(ts1);
        while(fcloverMap_it!=fcloverMap.end()&&fcloverMap_it->first<ts2){
            corrts = (Long64_t) fcloverMap_it->first;
            correntry = fcloverMap_it->second;
            if (corrts!=check_time){
                ftrGamma->GetEvent(correntry);
                check_time=corrts;
                CloverHit* hit=new CloverHit;
                fclover->Copy(*hit);
                flocalimp->AddGammaHit(hit);
                ncorr++;
            }
            fcloverMap_it++;
        }
        if (ncorr>0) ncorrwGamma++;

        /*
        //! correlated event with anc
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;

        fancMap_it = fancMap.lower_bound(ts1);
        while(fancMap_it!=fancMap.end()&&fancMap_it->first<ts2){
            corrts = (Long64_t) fancMap_it->first;
            correntry = fancMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalimp->AddAncHit(hit);
                ncorr++;
            }
            fancMap_it++;
        }
        if (ncorr>0) ncorrwAnc++;
        */

        //! correlated event with f11l
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapL_it = fF11LMap.lower_bound(ts1);
        while(fF11MapL_it!=fF11LMap.end()&&fF11MapL_it->first<ts2){
            corrts = (Long64_t) fF11MapL_it->first;
            correntry = fF11MapL_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalimp->GetF11Beam()->SetTF11L(fanc->GetTimestamp());
                flocalimp->GetF11Beam()->SetEF11L(fanc->GetEnergy());
                flocalimp->GetF11Beam()->SetNF11L(1);
                ncorr++;
                break;
            //}
            fF11MapL_it++;
        }
        if (ncorr>0) ncorrwF11L++;

        //! correlated event with f11r
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            correntry = fF11MapR_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalimp->GetF11Beam()->SetTF11R(fanc->GetTimestamp());
                flocalimp->GetF11Beam()->SetEF11R(fanc->GetEnergy());
                flocalimp->GetF11Beam()->SetNF11R(1);
                ncorr++;
                break;
            //}
            fF11MapR_it++;
        }
        if (ncorr>0) ncorrwF11R++;


        //! correlate event with vetotop
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoTopMap_it = fVetoTopMap.lower_bound(ts1);
        while(fVetoTopMap_it!=fVetoTopMap.end()&&fVetoTopMap_it->first<ts2){
            corrts = (Long64_t) fVetoTopMap_it->first;
            correntry = fVetoTopMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalimp->GetF11Beam()->SetTVetoTop(fanc->GetTimestamp());
                flocalimp->GetF11Beam()->SetEVetoTop(fanc->GetEnergy());
                flocalimp->GetF11Beam()->SetNVetoTop(1);
                ncorr++;
                break;
            }
            fVetoTopMap_it++;
        }
        if (ncorr>0) ncorrwVetoTop++;

        //! correlate event with vetobot
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoBotMap_it = fVetoBotMap.lower_bound(ts1);
        while(fVetoBotMap_it!=fVetoBotMap.end()&&fVetoBotMap_it->first<ts2){
            corrts = (Long64_t) fVetoBotMap_it->first;
            correntry = fVetoBotMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalimp->GetF11Beam()->SetTVetoBot(fanc->GetTimestamp());
                flocalimp->GetF11Beam()->SetEVetoBot(fanc->GetEnergy());
                flocalimp->GetF11Beam()->SetNVetoBot(1);
                ncorr++;
                break;
            }
            fVetoBotMap_it++;
        }
        if (ncorr>0) ncorrwVetoBot++;

        //! correlate event with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            correntry = fVetoDownMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalimp->GetF11Beam()->SetTVetoDown(fanc->GetTimestamp());
                flocalimp->GetF11Beam()->SetEVetoDown(fanc->GetEnergy());
                flocalimp->GetF11Beam()->SetNVetoDown(1);
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }
        if (ncorr>0) ncorrwVetoDown++;

        ftreeImplant->Fill();
    }
    cout<<"finished merging Implantation" <<endl;
    cout<<ncorrwBR<<" "<<ncorrwNeutron<<" "<<ncorrwGamma<<endl;
    cout<<ncorrwF11L<<" "<<ncorrwF11R<<" "<<ncorrwVetoTop<<" "<<ncorrwVetoBot<<" "<<ncorrwVetoDown<<endl;
}

void Merger::DoMergeBeta()
{
    unsigned int ncorrwGamma = 0;
    unsigned int ncorrwAnc = 0;
    unsigned int ncorrwF11R = 0;
    unsigned int ncorrwF11L = 0;
    unsigned int ncorrwVetoTop = 0;
    unsigned int ncorrwVetoBot = 0;
    unsigned int ncorrwVetoDown = 0;

    for (faidaBetaMap_it=faidaBetaMap.begin();faidaBetaMap_it!=faidaBetaMap.end();faidaBetaMap_it++){
        unsigned long long ts  = faidaBetaMap_it->first;
        unsigned int entry = faidaBetaMap_it->second.first;
        unsigned short clusterId = faidaBetaMap_it->second.second;
        ftrAIDABeta->GetEvent(entry);
        flocalbeta->Clear();

        //! fill beta clusters first
        faidaBeta->GetCluster(clusterId)->Copy(*flocalbeta->GetBeta());
        flocalbeta->SetTimeStamp(ts);
        //! and some commond stuff for light particles rejection
        flocalbeta->SetMult(faidaBeta->GetMult());
        flocalbeta->SetHitMultZ(faidaBeta->GetZHitMult());
        flocalbeta->SetNHitZ(faidaBeta->GetNHitZ());
        flocalbeta->SetMultX(faidaBeta->GetMultXs());
        flocalbeta->SetMultY(faidaBeta->GetMultYs());
        flocalbeta->SetNClusters(faidaBeta->GetNClusters());
        flocalbeta->SetClusterMultZ(faidaBeta->GetClustersMultZ());
        flocalbeta->SetNClusterZ(faidaBeta->GetNClustersZ());
        flocalbeta->SetMaxZ(faidaBeta->GetMaxZ());

        //! correlate event with gamma
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fBetaGammaTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fBetaGammaTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        fcloverMap_it = fcloverMap.lower_bound(ts1);
        while(fcloverMap_it!=fcloverMap.end()&&fcloverMap_it->first<ts2){
            corrts = (Long64_t) fcloverMap_it->first;
            correntry = fcloverMap_it->second;
            if (corrts!=check_time){
                ftrGamma->GetEvent(correntry);
                check_time=corrts;
                CloverHit* hit=new CloverHit;
                fclover->Copy(*hit);
                flocalbeta->AddGammaHit(hit);
                ncorr++;
            }
            fcloverMap_it++;
        }
        if (ncorr>0) ncorrwGamma++;

        /*
        //! correlated event with anc
        ts1 = (Long64_t)ts - (Long64_t)fBetaAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fancMap_it = fancMap.lower_bound(ts1);
        while(fancMap_it!=fancMap.end()&&fancMap_it->first<ts2){
            corrts = (Long64_t) fancMap_it->first;
            correntry = fancMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalbeta->AddAncHit(hit);
                ncorr++;
            }
            fancMap_it++;
        }
        if (ncorr>0) ncorrwAnc++;
        */


        //! correlated event with f11l
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapL_it = fF11LMap.lower_bound(ts1);
        while(fF11MapL_it!=fF11LMap.end()&&fF11MapL_it->first<ts2){
            corrts = (Long64_t) fF11MapL_it->first;
            correntry = fF11MapL_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalbeta->GetF11Beam()->SetTF11L(fanc->GetTimestamp());
                flocalbeta->GetF11Beam()->SetEF11L(fanc->GetEnergy());
                flocalbeta->GetF11Beam()->SetNF11L(1);
                ncorr++;
                break;
            //}
            fF11MapL_it++;
        }
        if (ncorr>0) ncorrwF11L++;

        //! correlated event with f11r
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            correntry = fF11MapR_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalbeta->GetF11Beam()->SetTF11R(fanc->GetTimestamp());
                flocalbeta->GetF11Beam()->SetEF11R(fanc->GetEnergy());
                flocalbeta->GetF11Beam()->SetNF11R(1);
                ncorr++;
                break;
            //}
            fF11MapR_it++;
        }
        if (ncorr>0) ncorrwF11R++;


        //! correlate event with vetotop
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoTopMap_it = fVetoTopMap.lower_bound(ts1);
        while(fVetoTopMap_it!=fVetoTopMap.end()&&fVetoTopMap_it->first<ts2){
            corrts = (Long64_t) fVetoTopMap_it->first;
            correntry = fVetoTopMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalbeta->GetF11Beam()->SetTVetoTop(fanc->GetTimestamp());
                flocalbeta->GetF11Beam()->SetEVetoTop(fanc->GetEnergy());
                flocalbeta->GetF11Beam()->SetNVetoTop(1);
                ncorr++;
                break;
            }
            fVetoTopMap_it++;
        }
        if (ncorr>0) ncorrwVetoTop++;

        //! correlate event with vetobot
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoBotMap_it = fVetoBotMap.lower_bound(ts1);
        while(fVetoBotMap_it!=fVetoBotMap.end()&&fVetoBotMap_it->first<ts2){
            corrts = (Long64_t) fVetoBotMap_it->first;
            correntry = fVetoBotMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalbeta->GetF11Beam()->SetTVetoBot(fanc->GetTimestamp());
                flocalbeta->GetF11Beam()->SetEVetoBot(fanc->GetEnergy());
                flocalbeta->GetF11Beam()->SetNVetoBot(1);
                ncorr++;
                break;
            }
            fVetoBotMap_it++;
        }
        if (ncorr>0) ncorrwVetoBot++;

        //! correlate event with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            correntry = fVetoDownMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalbeta->GetF11Beam()->SetTVetoDown(fanc->GetTimestamp());
                flocalbeta->GetF11Beam()->SetEVetoDown(fanc->GetEnergy());
                flocalbeta->GetF11Beam()->SetNVetoDown(1);
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }
        if (ncorr>0) ncorrwVetoDown++;

        ftreeBeta->Fill();
    }

    cout<<"finished merging Beta" <<endl;
    cout<<ncorrwGamma<<endl;
    cout<<ncorrwF11L<<" "<<ncorrwF11R<<" "<<ncorrwVetoTop<<" "<<ncorrwVetoBot<<" "<<ncorrwVetoDown<<endl;

}

void Merger::DoMergeNeutron()
{
    unsigned int ncorrwGamma = 0;
    unsigned int ncorrwAnc = 0;
    unsigned int ncorrwF11R = 0;
    unsigned int ncorrwF11L = 0;
    unsigned int ncorrwVetoTop = 0;
    unsigned int ncorrwVetoBot = 0;
    unsigned int ncorrwVetoDown = 0;

    for (fhe3Map_it=fhe3Map.begin();fhe3Map_it!=fhe3Map.end();fhe3Map_it++){
        unsigned long long ts  = fhe3Map_it->first;
        unsigned int entry = fhe3Map_it->second;
        ftrNeutron->GetEvent(entry);
        flocalneutron->Clear();

        //! fill neutron first
        fneutron->Copy(*flocalneutron->GetNeutron());
        flocalneutron->SetTimeStamp(ts);

        //! correlate event with gamma
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fNeuGammaTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fNeuGammaTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        fcloverMap_it = fcloverMap.lower_bound(ts1);
        while(fcloverMap_it!=fcloverMap.end()&&fcloverMap_it->first<ts2){
            corrts = (Long64_t) fcloverMap_it->first;
            correntry = fcloverMap_it->second;
            if (corrts!=check_time){
                ftrGamma->GetEvent(correntry);
                check_time=corrts;
                CloverHit* hit=new CloverHit;
                fclover->Copy(*hit);
                flocalneutron->AddGammaHit(hit);
                ncorr++;
            }
            fcloverMap_it++;
        }
        if (ncorr>0) ncorrwGamma++;

        /*
        //! correlate event with anc
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fancMap_it = fancMap.lower_bound(ts1);
        while(fancMap_it!=fancMap.end()&&fancMap_it->first<ts2){
            corrts = (Long64_t) fancMap_it->first;
            correntry = fancMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalneutron->AddAncHit(hit);
                ncorr++;
            }
            fancMap_it++;
        }
        if (ncorr>0) ncorrwAnc++;
        */

        //! correlated event with f11l
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapL_it = fF11LMap.lower_bound(ts1);
        while(fF11MapL_it!=fF11LMap.end()&&fF11MapL_it->first<ts2){
            corrts = (Long64_t) fF11MapL_it->first;
            correntry = fF11MapL_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalneutron->GetF11Beam()->SetTF11L(fanc->GetTimestamp());
                flocalneutron->GetF11Beam()->SetEF11L(fanc->GetEnergy());
                flocalneutron->GetF11Beam()->SetNF11L(1);
                ncorr++;
                break;
            //}
            fF11MapL_it++;
        }
        if (ncorr>0) ncorrwF11L++;

        //! correlated event with f11r
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            correntry = fF11MapR_it->second;
            //if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalneutron->GetF11Beam()->SetTF11R(fanc->GetTimestamp());
                flocalneutron->GetF11Beam()->SetEF11R(fanc->GetEnergy());
                flocalneutron->GetF11Beam()->SetNF11R(1);
                ncorr++;
                break;
            //}
            fF11MapR_it++;
        }
        if (ncorr>0) ncorrwF11R++;


        //! correlate event with vetotop
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoTopMap_it = fVetoTopMap.lower_bound(ts1);
        while(fVetoTopMap_it!=fVetoTopMap.end()&&fVetoTopMap_it->first<ts2){
            corrts = (Long64_t) fVetoTopMap_it->first;
            correntry = fVetoTopMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalneutron->GetF11Beam()->SetTVetoTop(fanc->GetTimestamp());
                flocalneutron->GetF11Beam()->SetEVetoTop(fanc->GetEnergy());
                flocalneutron->GetF11Beam()->SetNVetoTop(1);
                ncorr++;
                break;
            }
            fVetoTopMap_it++;
        }
        if (ncorr>0) ncorrwVetoTop++;

        //! correlate event with vetobot
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoBotMap_it = fVetoBotMap.lower_bound(ts1);
        while(fVetoBotMap_it!=fVetoBotMap.end()&&fVetoBotMap_it->first<ts2){
            corrts = (Long64_t) fVetoBotMap_it->first;
            correntry = fVetoBotMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalneutron->GetF11Beam()->SetTVetoBot(fanc->GetTimestamp());
                flocalneutron->GetF11Beam()->SetEVetoBot(fanc->GetEnergy());
                flocalneutron->GetF11Beam()->SetNVetoBot(1);
                ncorr++;
                break;
            }
            fVetoBotMap_it++;
        }
        if (ncorr>0) ncorrwVetoBot++;

        //! correlate event with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            correntry = fVetoDownMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalneutron->GetF11Beam()->SetTVetoDown(fanc->GetTimestamp());
                flocalneutron->GetF11Beam()->SetEVetoDown(fanc->GetEnergy());
                flocalneutron->GetF11Beam()->SetNVetoDown(1);
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }
        if (ncorr>0) ncorrwVetoDown++;

        ftreeNeutron->Fill();
    }
    cout<<"finished merging Neutrons" <<endl;
    cout<<ncorrwGamma<<endl;
    cout<<ncorrwF11L<<" "<<ncorrwF11R<<" "<<ncorrwVetoTop<<" "<<ncorrwVetoBot<<" "<<ncorrwVetoDown<<endl;
}
void Merger::DoMergeAnc()
{
    unsigned int ncorrwF11R = 0;
    unsigned int ncorrwVetoTop = 0;
    unsigned int ncorrwVetoBot = 0;
    unsigned int ncorrwVetoDown = 0;
    unsigned int ncorrwNeutron = 0;

    for (fF11MapL_it=fF11LMap.begin();fF11MapL_it!=fF11LMap.end();fF11MapL_it++){
        unsigned long long ts  = fF11MapL_it->first;
        unsigned int entry = fF11MapL_it->second;
        ftrAnc->GetEvent(entry);
        flocalancillary->Clear();
        flocalancillary->SetEF11L(fanc->GetEnergy());
        flocalancillary->SetTF11L(fanc->GetTimestamp());
        flocalancillary->SetNF11L(1);

        //! correlate event with F11Rs
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fF11LRTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fF11LRTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            correntry = fF11MapR_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalancillary->SetTF11R(fanc->GetTimestamp());
                flocalancillary->SetEF11R(fanc->GetEnergy());
                flocalancillary->SetNF11R(1);
                ncorr++;
                break;
            }
            fF11MapR_it++;
        }
        if (ncorr>0) ncorrwF11R++;

        //! correlate event with vetotop
        ts1 = (Long64_t)ts - (Long64_t)fF11LVetoTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fF11LVetoTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoTopMap_it = fVetoTopMap.lower_bound(ts1);
        while(fVetoTopMap_it!=fVetoTopMap.end()&&fVetoTopMap_it->first<ts2){
            corrts = (Long64_t) fVetoTopMap_it->first;
            correntry = fVetoTopMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalancillary->SetTVetoTop(fanc->GetTimestamp());
                flocalancillary->SetEVetoTop(fanc->GetEnergy());
                flocalancillary->SetNVetoTop(1);
                ncorr++;
                break;
            }
            fVetoTopMap_it++;
        }
        if (ncorr>0) ncorrwVetoTop++;

        //! correlate event with vetobot
        ts1 = (Long64_t)ts - (Long64_t)fF11LVetoTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fF11LVetoTWup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoBotMap_it = fVetoBotMap.lower_bound(ts1);
        while(fVetoBotMap_it!=fVetoBotMap.end()&&fVetoBotMap_it->first<ts2){
            corrts = (Long64_t) fVetoBotMap_it->first;
            correntry = fVetoBotMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalancillary->SetTVetoBot(fanc->GetTimestamp());
                flocalancillary->SetEVetoBot(fanc->GetEnergy());
                flocalancillary->SetNVetoBot(1);
                ncorr++;
                break;
            }
            fVetoBotMap_it++;
        }
        if (ncorr>0) ncorrwVetoBot++;

        //! correlate event with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fF11LVetoDownlow;
        ts2 = (Long64_t)ts + (Long64_t)fF11LVetoDownup;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            correntry = fVetoDownMap_it->second;
            if (corrts!=check_time){
                ftrAnc->GetEvent(correntry);
                check_time=corrts;
                flocalancillary->SetTVetoDown(fanc->GetTimestamp());
                flocalancillary->SetEVetoDown(fanc->GetEnergy());
                flocalancillary->SetNVetoDown(1);
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }
        if (ncorr>0) ncorrwVetoDown++;

        //! correlated event with neutron
        ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWup;
        ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWlow;
        ncorr = 0;
        corrts = 0;
        correntry = 0;
        check_time =0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            correntry = fhe3Map_it->second;
            if (corrts!=check_time){
                ftrNeutron->GetEvent(correntry);
                check_time=corrts;
                BELENHit* hit=new BELENHit;
                fneutron->Copy(*hit);
                flocalancillary->AddNeutronHit(hit);
                ncorr++;
            }
            fhe3Map_it++;
        }
        if (ncorr>0) ncorrwNeutron++;

        ftreeAncillary->Fill();
    }
    cout<<"Finished merging F11L" <<endl;
    cout<<ncorrwF11R<<" "<<ncorrwVetoTop<<" "<<ncorrwVetoBot<<" "<<ncorrwVetoDown<<" "<<ncorrwNeutron<<endl;
}
