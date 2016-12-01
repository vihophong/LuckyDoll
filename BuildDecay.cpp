#include "BuildDecay.h"
#include "TROOT.h"
BuildDecay::BuildDecay()
{
    fMode = 0;
    fBetaImplantTWup = 1e9;
    fBetaImplantTWlow = 6e9;
    fBetaNeutronTWup = 5e5;
    fBetaNeutronTWlow = 1e5;
    fmaxdeltaxy = 2.;

    fBetaF11BeamTWup = 40000;
    fBetaF11BeamTWlow = 50000;
    fIonF11BeamTWup = 40000;
    fIonF11BeamTWlow = 50000;
}

BuildDecay::~BuildDecay()
{
    delete fh1;
}

void BuildDecay::Init()
{
    ftreedecay = 0;
    fimplant = NULL;
    fbeta = NULL;
    fneutron = NULL;
    fnF11Beam = 0;
    fF11Beam = NULL;

    flocalneutron = new Neutrons;

    fimplanti = 0;
    fnentriesBeta = 0;
    fnentriesImp = 0;
    fnentriesNeutron = 0;
    //! init input files

    finputFile = new TFile(finput);
    finputFile->GetObject("implant",ftrImp);
    finputFile->GetObject("beta",ftrBeta);
    finputFile->GetObject("neutron",ftrNeutron);
    finputFile->GetObject("beam",ftrF11Beam);
    ftrImp->SetBranchAddress("implant",&fimplant);
    ftrBeta->SetBranchAddress("beta",&fbeta);
    ftrNeutron->SetBranchAddress("neutron",&fneutron);
    ftrF11Beam->SetBranchAddress("beam",&fF11Beam);

    flocalbeta = fbeta;
    flocalimp = fimplant;
    flocalF11Beam = new Ancillary;

    fnentriesImp = ftrImp->GetEntries();
    fnentriesBeta = ftrBeta->GetEntries();
    fnentriesNeutron = ftrNeutron->GetEntries();
    fnentriesF11Beam = ftrF11Beam->GetEntries();

    cout<<"Reading "<<fnentriesImp<<" implants, "<<fnentriesBeta<<" betas and "<<fnentriesNeutron<<" neutrons "<<fnentriesF11Beam<<" f11beam"<<endl;
    cout<<"Printing first few timestamps:"<<endl;
    for (unsigned int i=0;i<10;i++){
        ftrImp->GetEvent(i);
        ftrBeta->GetEvent(i);
        ftrNeutron->GetEvent(i);
        ftrF11Beam->GetEvent(i);
        cout<<"implantts "<<fimplant->GetTimeStamp()<<" - betats "
           <<fbeta->GetTimeStamp()<<" -neutronts "<<fneutron->GetTimeStamp()
          <<" -f11beamts "<<fF11Beam->GetTimeStamp()<<endl;
    }

    fh1=new TH1F("h1","h1",500,0,100e3);
}

void BuildDecay::BookTree(TTree* treeDecay){
    ftreedecay = treeDecay;
    if (fMode==1){
        ftreedecay->Branch("fnF11Beam",&fnF11Beam,"fnF11Beam/s");
        ftreedecay->Branch("beta",&flocalbeta);
        ftreedecay->Branch("neutron",&flocalneutron);
        ftreedecay->Branch("f11beam",&flocalF11Beam);
    }else if (fMode==2){
        ftreedecay->Branch("fnF11Beam",&fnF11Beam,"fnF11Beam/s");
        ftreedecay->Branch("implant",&flocalimp);
        ftreedecay->Branch("neutron",&flocalneutron);
        ftreedecay->Branch("f11beam",&flocalF11Beam);
    }else{
        ftreedecay->Branch("fdeltaxy",&fdeltaxy,"fdeltaxy/D");
        ftreedecay->Branch("fimplanti",&fimplanti,"fimplanti/s");
        ftreedecay->Branch("implant",&flocalimp);
        ftreedecay->Branch("beta",&flocalbeta);
        ftreedecay->Branch("neutron",&flocalneutron);
    }
    ftreedecay->BranchRef();
}

void BuildDecay::ReadImplant()
{
    fimplantMap.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesImp; jentry++){
        ftrImp->GetEvent(jentry);
        //! make gates and selections here
        aidaSimpleStruct aida;
        aida.x = fimplant->GetIon()->GetHitPositionX();
        aida.y = fimplant->GetIon()->GetHitPositionY();
        aida.z = (unsigned short)fimplant->GetIon()->GetHitPositionZ();
        fimplantMap.insert(make_pair(fimplant->GetTimeStamp(),make_pair(jentry,aida)));
    }
    cout<<"Finished reading implantation table with "<<fimplantMap.size()<<" rows"<<endl;
}
void BuildDecay::ReadImplant2()
{
    fimplantMap2.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesImp; jentry++){
        ftrImp->GetEvent(jentry);
        fimplantMap2.insert(make_pair(fimplant->GetTimeStamp(),jentry));
    }
    cout<<"Finished reading implantation table with "<<fimplantMap2.size()<<" rows"<<endl;
}

void BuildDecay::ReadBeta()
{
    fbetaMap.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesBeta; jentry++){
        ftrBeta->GetEvent(jentry);
        //! make gates and selections here
        fbetaMap.insert(make_pair(fbeta->GetTimeStamp(),jentry));
    }
    cout<<"Finished reading beta table with "<<fbetaMap.size()<<" rows"<<endl;
}
void BuildDecay::ReadBeta2()
{
    fbetaMap2.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesBeta; jentry++){
        ftrBeta->GetEvent(jentry);
        //! make gates and selections here
        aidaSimpleStruct aida;
        aida.x = fbeta->GetBeta()->GetHitPositionX();
        aida.y = fbeta->GetBeta()->GetHitPositionY();
        aida.z = (unsigned short)fbeta->GetBeta()->GetHitPositionZ();
        fbetaMap2.insert(make_pair(fbeta->GetTimeStamp(),make_pair(jentry,aida)));
    }
    cout<<"Finished reading beta table with "<<fbetaMap2.size()<<" rows"<<endl;
}

void BuildDecay::ReadNeutron()
{
    fneutronMap.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesNeutron; jentry++){
        ftrNeutron->GetEvent(jentry);
        //! make gates and selections here
        fneutronMap.insert(make_pair(fneutron->GetTimeStamp(),jentry));
    }
    cout<<"Finished reading neutron table with "<<fneutronMap.size()<<" rows"<<endl;
}

void BuildDecay::ReadF11Beam()
{
    fF11BeamMap.clear();
    for (unsigned int jentry = 0; jentry< (unsigned int) fnentriesF11Beam; jentry++){
        ftrF11Beam->GetEvent(jentry);
        //! make gates and selections here
        fF11BeamMap.insert(make_pair(fF11Beam->GetTimeStamp(),jentry));
    }
    cout<<"Finished reading F11beam table with "<<fF11BeamMap.size()<<" rows"<<endl;
}


void BuildDecay::DoBuildDecay()
{
    cout<<"\n\n***********Start Building Decay************ \n\n"<<endl;
    unsigned int ncorrwIon = 0;
    unsigned int ncorrwNeutron = 0;
    flocalneutron->Clear();
    fimplanti = 0;
    fdeltaxy = -1;
    Int_t k = 0;
    for (fbetaMap_it=fbetaMap.begin();fbetaMap_it!=fbetaMap.end();fbetaMap_it++){
        if (k%5000==0) cout<<(Double_t)k/(Double_t)fnentriesBeta*100.<<" % completed, file size ="<<(Double_t)foutfile->GetSize()/(Double_t)1000000.<<" MB \r"<<flush;
        fimplanti = 0;
        fdeltaxy = -1;       
        unsigned long long ts  = fbetaMap_it->first;
        unsigned int entry = fbetaMap_it->second;
        ftrBeta->GetEvent(entry);
        double betax = fbeta->GetBeta()->GetHitPositionX();
        double betay = fbeta->GetBeta()->GetHitPositionY();
        unsigned short betaz = (unsigned short) fbeta->GetBeta()->GetHitPositionZ();
        //! Clear neutron
        flocalneutron->Clear();
        //! correlate event with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fBetaNeutronTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fBetaNeutronTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;



        fneutronMap_it = fneutronMap.lower_bound(ts1);
        while(fneutronMap_it!=fneutronMap.end()&&fneutronMap_it->first<ts2){
            corrts = (Long64_t) fneutronMap_it->first;
            correntry = fneutronMap_it->second;
            if (corrts!=check_time){
                ftrNeutron->GetEvent(correntry);
                check_time=corrts;
                Neutron* hit=new Neutron;
                fneutron->Copy(*hit);
                flocalneutron->AddNeutron(hit);
                ncorr++;
            }
            fneutronMap_it++;
        }
        if (ncorr>0) ncorrwNeutron++;

        //! correlate event with implant
        ts1 = (Long64_t)ts - (Long64_t)fBetaImplantTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaImplantTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        fimplantMap_it = fimplantMap.lower_bound(ts1);
        while(fimplantMap_it!=fimplantMap.end()&&fimplantMap_it->first<ts2){
            corrts = (Long64_t) fimplantMap_it->first;
            correntry = fimplantMap_it->second.first;
            aidaSimpleStruct ion = fimplantMap_it->second.second;
            if (corrts!=check_time&&ion.z==betaz&&ion.x>=0&&ion.y>=0){
                fdeltaxy = sqrt((ion.x-betax)*(ion.x-betax)+(ion.y-betay)*(ion.y-betay));
                if (fdeltaxy>fmaxdeltaxy){
                    fimplantMap_it++;
                    continue;
                }
                ftrImp->GetEntry(correntry);
                //! fill data here!
                fimplanti ++;
                ftreedecay->Fill();
            }
            fimplantMap_it++;
        }
        if (fimplanti) ncorrwIon++;
        //ftreedecay->Fill();
        k++;

    }
    cout<<"finished build Decay!"<<endl;
    cout<<ncorrwNeutron<<" "<<ncorrwIon<<endl;
}

void BuildDecay::DoBuildDecay3()
{
    cout<<"\n\n***********Start Building Beta-Neutron************ \n\n"<<endl;
    unsigned int ncorrwNeutron = 0;
    unsigned int ncorrwF11Beam = 0;
    flocalneutron->Clear();
    Int_t k = 0;
    for (fbetaMap_it=fbetaMap.begin();fbetaMap_it!=fbetaMap.end();fbetaMap_it++){
        if (k%5000==0) cout<<(Double_t)k/(Double_t)fnentriesBeta*100.<<" % completed, file size ="<<(Double_t)foutfile->GetSize()/(Double_t)1000000.<<" MB \r"<<flush;
        unsigned long long ts  = fbetaMap_it->first;
        unsigned int entry = fbetaMap_it->second;
        ftrBeta->GetEvent(entry);
        //! Clear neutron
        flocalneutron->Clear();
        fnF11Beam = 0;
        flocalF11Beam->Clear();
        //! correlate event with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fBetaNeutronTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fBetaNeutronTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fneutronMap_it = fneutronMap.lower_bound(ts1);
        while(fneutronMap_it!=fneutronMap.end()&&fneutronMap_it->first<ts2){
            corrts = (Long64_t) fneutronMap_it->first;
            correntry = fneutronMap_it->second;
            if (corrts!=check_time){
                ftrNeutron->GetEvent(correntry);
                check_time=corrts;
                Neutron* hit=new Neutron;
                fneutron->Copy(*hit);
                flocalneutron->AddNeutron(hit);
                ncorr++;
            }
            fneutronMap_it++;
        }
        if (ncorr>0) ncorrwNeutron++;

        //! correlate event with implant
        ts1 = (Long64_t)ts - (Long64_t)fBetaF11BeamTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaF11BeamTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;

        fF11BeamMap_it = fF11BeamMap.lower_bound(ts1);
        while(fF11BeamMap_it!=fF11BeamMap.end()&&fF11BeamMap_it->first<ts2){
            corrts = (Long64_t) fF11BeamMap_it->first;
            correntry = fF11BeamMap_it->second;
            if (corrts!=check_time){
                ftrF11Beam->GetEvent(correntry);
                fF11Beam->Copy(*flocalF11Beam);
                fnF11Beam++;
                check_time=corrts;
                ncorr++;
                break;
            }
            fF11BeamMap_it++;
        }
        if (ncorr>0) ncorrwF11Beam++;

        ftreedecay->Fill();
        k++;
    }
    cout<<"finished build beta-neutron!"<<endl;
    cout<<ncorrwNeutron<<" "<<ncorrwF11Beam<<endl;
}

void BuildDecay::DoBuildDecay4()
{
    cout<<"\n\n***********Start Building Ion-Neutron************ \n\n"<<endl;
    unsigned int ncorrwNeutron = 0;
    unsigned int ncorrwF11Beam = 0;
    flocalneutron->Clear();
    Int_t k = 0;
    for (fimplantMap_it=fimplantMap.begin();fimplantMap_it!=fimplantMap.end();fimplantMap_it++){
        if (k%5000==0) cout<<(Double_t)k/(Double_t)fnentriesImp*100.<<" % completed, file size ="<<(Double_t)foutfile->GetSize()/(Double_t)1000000.<<" MB \r"<<flush;
        unsigned long long ts  = fimplantMap_it->first;
        unsigned int entry = fimplantMap_it->second.first;
        ftrImp->GetEvent(entry);
        //! Clear neutron
        flocalneutron->Clear();
        fnF11Beam = 0;
        flocalF11Beam->Clear();
        //! correlate event with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fBetaNeutronTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fBetaNeutronTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fneutronMap_it = fneutronMap.lower_bound(ts1);
        while(fneutronMap_it!=fneutronMap.end()&&fneutronMap_it->first<ts2){
            corrts = (Long64_t) fneutronMap_it->first;
            correntry = fneutronMap_it->second;
            if (corrts!=check_time){
                ftrNeutron->GetEvent(correntry);
                check_time=corrts;
                Neutron* hit=new Neutron;
                fneutron->Copy(*hit);
                flocalneutron->AddNeutron(hit);
                ncorr++;
            }
            fneutronMap_it++;
        }
        if (ncorr>0) ncorrwNeutron++;

        //! correlate event with implant
        ts1 = (Long64_t)ts - (Long64_t)fBetaF11BeamTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaF11BeamTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;

        fF11BeamMap_it = fF11BeamMap.lower_bound(ts1);
        while(fF11BeamMap_it!=fF11BeamMap.end()&&fF11BeamMap_it->first<ts2){
            corrts = (Long64_t) fF11BeamMap_it->first;
            correntry = fF11BeamMap_it->second;
            if (corrts!=check_time){
                ftrF11Beam->GetEvent(correntry);
                fF11Beam->Copy(*flocalF11Beam);
                fnF11Beam++;
                check_time=corrts;
                ncorr++;
                break;
            }
            fF11BeamMap_it++;
        }
        if (ncorr>0) ncorrwF11Beam++;

        ftreedecay->Fill();
        k++;
    }
    cout<<"finished build ion-neutron!"<<endl;
    cout<<ncorrwNeutron<<" "<<ncorrwF11Beam<<endl;
}

void BuildDecay::DoBuildDecay2()
{
    cout<<"\n\n***********Start Building Decay************ \n\n"<<endl;
    unsigned int ncorrwIon = 0;
    unsigned int ncorrwNeutron = 0;
    flocalneutron->Clear();
    fimplanti = 0;
    fdeltaxy = -1;
    Int_t k = 0;
    for (fimplantMap_it2=fimplantMap2.begin();fimplantMap_it2!=fimplantMap2.end();fimplantMap_it2++){
        if (k%10000==0) cout<<k<<"-"<<fnentriesImp<<" \r "<<flush;
        unsigned long long ts  = fimplantMap_it2->first;
        unsigned int entry = fimplantMap_it2->second;
        ftrImp->GetEvent(entry);
        fimplanti = 0;
        fdeltaxy = -1;
        double impx = fimplant->GetIon()->GetHitPositionX();
        double impy = fimplant->GetIon()->GetHitPositionY();
        unsigned short impz = (unsigned short) fimplant->GetIon()->GetHitPositionZ();

        //! correlate event with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fBetaImplantTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fBetaImplantTWup;
        Long64_t ncorr = 0;
        Long64_t corrts = 0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fbetaMap_it2 = fbetaMap2.lower_bound(ts1);
        while(fbetaMap_it2!=fbetaMap2.end()&&fbetaMap_it2->first<ts2){
            corrts = (Long64_t) fbetaMap_it2->first;
            correntry = fbetaMap_it2->second.first;
            aidaSimpleStruct beta = fbetaMap_it2->second.second;
            if (corrts!=check_time&&beta.z==impz&&beta.x>=0&&beta.y>=0){
                fdeltaxy = sqrt((beta.x-impx)*(beta.x-impx)+(beta.y-impy)*(beta.y-impy));
                if (fdeltaxy>fmaxdeltaxy){
                    fbetaMap_it2++;
                    continue;
                }
                ftrBeta->GetEntry(correntry);
                //! fill data here!
                fimplanti ++;
                ftreedecay->Fill();
            }
            fbetaMap_it2++;
        }
        if (fimplanti) ncorrwIon++;
        k++;
    }

    cout<<"finished build Decay!"<<endl;
    cout<<ncorrwNeutron<<" "<<ncorrwIon<<endl;
}
