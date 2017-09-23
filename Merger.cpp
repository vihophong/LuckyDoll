#include "Merger.h"
Merger::Merger():fbigrips()
{
    fmaxmult=400;
    fmaxnpixels=3.;
    fIonBetaTWlow=10000000000;
    fIonBetaTWup=10000000000;
    fIonPidTWup = 0;
    fIonPidTWlow = 50000;
    fIonNeutronTWup = 800000;
    fIonNeutronTWlow = 100000;
    fIonGammaTWup = 5000;
    fIonGammaTWlow = 70000;
    fIonAncTWup = 10000;
    fIonAncTWlow = 80000;

    fIondETWup = 10000;
    fIondETWlow = 80000;

    fNeuGammaTWup = 100000;
    fNeuGammaTWlow = 800000;

    fNeuAncTWup = 600000;
    fNeuAncTWlow = 0;

    //fNeuAncTWup = 1000000;
    //fNeuAncTWlow = 1000000;

    fNeuBetaTWup = 600000;
    fNeuBetaTWlow = 600000;
    fNeuBetaoffset = -22000;
    //fNeuBetaoffset = 0;//for simulation
    fNeuBetaTWup= fNeuBetaTWup-fNeuBetaoffset;
    fNeuBetaTWlow= fNeuBetaTWlow+fNeuBetaoffset;//this is upper :)

    fBetaGammaTWup = 0;
    fBetaGammaTWlow = 30000;
    fBetaAncTWup = 10000;
    fBetaAncTWlow = 40000;

    fF11LRTWup = 400;
    fF11LRTWlow = 400;

    fF11LdEtopWup = 400;
    fF11LdEbotWup = 400;

    fF11LVetoTWup = 400;
    fF11LVetoTWlow  = 400;
    fF11LVetoDownup  = 3000;
    fF11LVetoDownlow  = 500;

    fF11LGammaTWlow = 1000;
    fF11LGammaTWup = 1000;

    finputAida=NULL;
    finputBigrips=NULL;
    finputBriken=NULL;

    fminneue = 100;
    fmaxneue = 1000;

    fmineneuvetodown=0.5;
    fmineneuf11=0;
}

Merger::~Merger()
{
    delete fh1;
}


void Merger::ReadPID(char* pidfile, Int_t ncutpts){
    std::ifstream ifspid(pidfile);
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;
    TString tempriname,tempria;
    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;
    }
     //! Setup PID cut
    for (Int_t i=0;i<nri;i++){
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);

        cutg[i]=new TCutG(nameri[i],ncutpts);
        cutg[i]->SetTitle(nameri[i]);
        cutg[i]->SetVarX("decay.aoq");
        cutg[i]->SetVarY("decay.zet");

        Double_t theta=0;
        Double_t x1=parmsri[i][0];
        Double_t y1=parmsri[i][1];
        Double_t r1=parmsri[i][2];
        Double_t r2=parmsri[i][3];

        Double_t phi1 = TMath::Min(0,360);
        Double_t phi2 = TMath::Max(0,360);
        Double_t kPI = 3.14159265358979323846;
        Double_t angle,dx,dy;

        Double_t x[ncutpts], y[ncutpts];
        Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
        Double_t ct   = TMath::Cos(kPI*theta/180);
        Double_t st   = TMath::Sin(kPI*theta/180);

        for (Int_t j=0;j<ncutpts;j++) {
           angle = phi1*kPI/180 + Double_t(j)*dphi;
           dx    = r1*TMath::Cos(angle);
           dy    = r2*TMath::Sin(angle);
           x[j]  = x1 + dx*ct - dy*st;
           y[j]  = y1 + dx*st + dy*ct;
           //cout<<x[j]<<"<"<<y[j]<<endl;
           cutg[i]->SetPoint(j,x[j],y[j]);
        }
        cutg[i]->SetPoint(ncutpts,x[0],y[0]);
        cutg[i]->SetName((char*)nameri[i].Data());
    }

    for (Int_t i=0;i<nri;i++){
        ftreeRI[i]=new TTree(Form("tree%s",(char*)nameri[i].Data()),Form("tree%s",(char*)nameri[i].Data()));
    }
}

void Merger::InitPIDSep()
{
    flocalimp=NULL;
    flocalbetaS = new IonBeta;
    flocalbetaS->Clear();
    //! init aida
    fmergedFile = new TFile(finputMerged);
    fmergedFile->GetObject("tree",ftrMerged);
    ftrMerged->SetBranchAddress("idx",&nionbetacorr);
    ftrMerged->SetBranchAddress("ion",&flocalimp);
    ftrMerged->GetBranch("beta")->SetAutoDelete(kFALSE);
    ftrMerged->SetBranchAddress("beta",&flocalbetaS);

    fnentriesMerged = ftrMerged->GetEntries();
    cout<<"Reading "<<fnentriesMerged<<" enetries in Merged tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrMerged->GetEvent(i);
        cout<<"ionts= "<<flocalimp->GetTimestamp()<<endl;
        cout<<"ionts= "<<flocalbetaS->GetTimestamp()<<endl;
    }

    TVectorD* deadtimecontainer=(TVectorD*) fmergedFile->Get(Form("deadtime"));
    Double_t* deadtimearray=deadtimecontainer->GetMatrixArray();
    fvetodeadtime=deadtimearray[0];
    fvetototaltime=deadtimearray[1];
    ff11vetodeadtime=deadtimearray[2];
    ff11vetototaltime=deadtimearray[3];
    fdownstreamvetodeadtime=deadtimearray[4];
    fdownstreamvetototaltime=deadtimearray[5];
    cout<<"F11R veto: deadtime="<<ff11vetodeadtime<<"\t---total time---\t"<<ff11vetototaltime<<endl;
    cout<<"Downstream veto:deadtime="<<fdownstreamvetodeadtime<<"\t---total time---\t"<<fdownstreamvetototaltime<<endl;
    cout<<"Final veto: deadtime="<<fvetodeadtime<<"\t---total time---\t"<<fvetototaltime<<endl;
}

void Merger::Init()
{
    ftree = 0;
    ftreeImplant = 0;
    ftreeBeta = 0;
    ftreeNeutron = 0;   


    flocalbeta = new TClonesArray("IonBeta");
    flocalbeta->Clear("C");
    flocalbetaS = new IonBeta;
    flocalbetaS->Clear();

    flocalimp = new IonBetaMult;
    flocalimp->Clear();
    nionbetacorr=0;

    faida = NULL;
    //faidaIon = NULL;
    //faidaBeta = NULL;
    fbigrips = NULL;
    fclover = NULL;
    fneutron = NULL;
    fanc = NULL;

    fnentriesAIDA=0;
    fnentriesAIDAIon = 0;
    fnentriesAIDABeta = 0;
    fnentriesBigrips = 0;
    fnentriesGamma = 0;
    fnentriesNeutron = 0;
    fnentriesAnc = 0;

    /*
    //! init aida
    if (finputAida!=NULL){
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
    }
    */

    //! init aida
    if (finputAida!=NULL){
        fAidaFile = new TFile(finputAida);
        fAidaFile->GetObject("aida",ftrAIDA);
        ftrAIDA->SetBranchAddress("aida",&faida);
        fnentriesAIDA = ftrAIDA->GetEntries();
        cout<<"Reading "<<fnentriesAIDA<<" enetries in AIDA tree"<<endl;
        cout<<"Printing first few timestamp:"<<endl;
        for (Long64_t i=0;i<10;i++){
            ftrAIDA->GetEvent(i);
            cout<<"ts= "<<faida->GetTimestamp()<<endl;
        }
    }


    if (finputBriken!=NULL){
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
    }


    if (finputBigrips!=NULL){
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
            cout<<"bigripsts "<<fbigrips->ts<<endl;
        }
        cout<<"\n\t\t\t****************************************\n"<<endl;
    }

    fh1=new TH1F("h1","h1",5000,-1000,1000);
}



void Merger::BookTreeTClone(TTree* tree,TTree* treemlh,TTree* treemlhp1n,TTree* treemlhp2n,TTree* treemlhp1nb,TTree* treemlhp2nb){
    ftree = tree;
    ftree->Branch("idx",&nionbetacorr,"idx/I");
    ftree->Branch("ion",&flocalimp);
    ftree->Branch("beta",&flocalbeta);
    ftree->BranchRef();
}

void Merger::BookTreeSingle(TTree* tree){
    ftree = tree;
    ftree->Branch("idx",&nionbetacorr,"idx/I");
    ftree->Branch("ion",&flocalimp);
    ftree->Branch("beta",&flocalbetaS);
    ftree->BranchRef();
}

void Merger::BookPIDSepTree(){
    ftree=new TTree("tree","tree");
    ftree->Branch("idx",&nionbetacorr,"idx/I");
    ftree->Branch("ion",&flocalimp);
    ftree->Branch("beta",&flocalbetaS);
    ftree->BranchRef();
    for (Int_t i=0;i<nri;i++){
        ftreeRI[i]->Branch("idx",&nionbetacorr,"idx/I");
        ftreeRI[i]->Branch("ion",&flocalimp);
        ftreeRI[i]->Branch("beta",&flocalbetaS);
        ftreeRI[i]->BranchRef();
    }
}

void Merger::ReadAIDA(unsigned int start, unsigned int stop)
{
    //! read ion
    unsigned int sstartI,sstopI;
    if (start==0) sstartI = 0; else sstartI = start;
    if (stop==0) sstopI = (unsigned int) fnentriesAIDA; else sstopI = stop;

    for (unsigned int jentry = sstartI;jentry < sstopI;jentry++){
        ftrAIDA->GetEvent(jentry);
        if (faida->GetID()==4){
            AIDASimpleStruct* aidasimple=new AIDASimpleStruct;
            faida->Copy(*aidasimple);
            faidaIonMap.insert(make_pair(faida->GetTimestamp(),aidasimple));
        }else if (faida->GetID()==5&&faida->GetMultiplicity()<fmaxmult){
            AIDASimpleStruct* aidasimple=new AIDASimpleStruct;
            faida->Copy(*aidasimple);
            faidaBetaMap.insert(make_pair(faida->GetTimestamp(),aidasimple));
        }
    }
    cout<<"Finished reading AIDA time table with "<<faidaIonMap.size()<<" entries for ion, and "<<faidaBetaMap.size()<<" entries for beta"<<endl;
}

void Merger::ReadBigrips()
{
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesBigrips;jentry++){
        ftrBigrips->GetEvent(jentry);
        fbigripsMap.insert(make_pair(fbigrips->ts,jentry));
    }
    cout<<"Finished reading bigrips ts table with "<<fbigripsMap.size()<<" rows"<<endl;
}

void Merger::ReadBRIKEN(unsigned int startN, unsigned int stopN,unsigned int startG, unsigned int stopG,unsigned int startA, unsigned int stopA)
{
    //! read neutron
    unsigned int sstartN,sstopN;
    if (startN==0) sstartN = 0; else sstartN = startN;
    if (stopN==0) sstopN = (unsigned int) fnentriesNeutron; else sstopN = stopN;

    for (unsigned int jentry = sstartN;jentry < sstopN;jentry++){
        ftrNeutron->GetEvent(jentry);
        BELENHit* hit=new BELENHit;
        fneutron->Copy(*hit);
        if (fneutron->GetEnergy()<fmaxneue&&fneutron->GetEnergy()>fminneue) fhe3Map.insert(make_pair(fneutron->GetTimestamp(),hit));
    }
    cout<<"Finished reading neutron  ts table with "<<fhe3Map.size()<<" rows"<<endl;

    //! read gamma
    unsigned int sstartG,sstopG;
    if (startG==0) sstartG = 0; else sstartG = startG;
    if (stopG==0) sstopG = (unsigned int) fnentriesGamma; else sstopG = stopG;

    for (unsigned int jentry = sstartG;jentry < sstopG;jentry++){
        ftrGamma->GetEvent(jentry);
        fcloverMap.insert(make_pair(fclover->GetTimestamp(),jentry));
    }
    cout<<"Finished reading gamma  ts table with "<<fcloverMap.size()<<" rows"<<endl;

    //! read anc
    unsigned int sstartA,sstopA;
    if (startA==0) sstartA = 0; else sstartA = startA;
    if (stopA==0) sstopA = (unsigned int) fnentriesAnc; else sstopA = stopA;

    for (unsigned int jentry = sstartA;jentry < sstopA;jentry++){
        ftrAnc->GetEvent(jentry);
        //fancMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==1) fF11RMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==2) fF11LMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==1) fVetoTopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==2) fVetoBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==3&&fanc->GetID()==1) fdETopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==3&&fanc->GetID()==2) fdEBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==4) fVetoDownMap.insert(make_pair(fanc->GetTimestamp(),jentry));
    }
    //cout<<"Finished reading anc  ts table with "<<fancMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11R  ts table with "<<fF11RMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11L  ts table with "<<fF11LMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Top ts table with "<<fVetoTopMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Bottom ts table with "<<fVetoBotMap.size()<<" rows"<<endl;
    cout<<"Finished reading dE Top ts table with "<<fdETopMap.size()<<" rows"<<endl;
    cout<<"Finished reading dE Bottom ts table with "<<fdEBotMap.size()<<" rows"<<endl;

    cout<<"Finished reading Veto Down ts table with "<<fVetoDownMap.size()<<" rows"<<endl;
}


void Merger::DoMergeTClone()
{
    Int_t k=0;
    Int_t ktotal=faidaIonMap.size();
    Int_t ncorrwbigrips=0;
    Int_t ncorrwbeta=0;

    for (faidaIonMap_it=faidaIonMap.begin();faidaIonMap_it!=faidaIonMap.end();faidaIonMap_it++){
        flocalimp->Clear();
        flocalbeta->Clear("C");
        //flocalbeta->Clear();
        if (k%1000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorrwbeta<<"\t ncorr w bigrips "<<ncorrwbigrips<<endl;
        //if (k>500) break;
        unsigned long long ts=faidaIonMap_it->first;
        AIDASimpleStruct* ion=(AIDASimpleStruct*) faidaIonMap_it->second;
        double impx=ion->GetHitPositionX();
        double impy=ion->GetHitPositionY();
        short impz=ion->GetHitPositionZ()+ion->GetDZ();//! correct dZ

        //! fill imp data here
        //AIDASimpleStruct* iontofill = flocalimp->GetIonBetaCluster();
        //ion->Copy(*iontofill);
        flocalimp->CopyFromAIDA(ion);
        //! Correlate imp with bigrips
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fbigripsMap_it = fbigripsMap.lower_bound(ts1);
        while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
            corrts = (Long64_t) fbigripsMap_it->first;
            correntry = fbigripsMap_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill data here
                ftrBigrips->GetEvent(correntry);
                flocalimp->AddBeam(*fbigrips);
                ncorr++;
                break;
            }
            fbigripsMap_it++;
        }
        if (ncorr>0) ncorrwbigrips++;

        //! Correlate with beta
        ts1 = (Long64_t)ts - fIonBetaTWlow;
        ts2 = (Long64_t)ts + fIonBetaTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        nionbetacorr=0;
        faidaBetaMap_it = faidaBetaMap.lower_bound(ts1);
        while(faidaBetaMap_it!=faidaBetaMap.end()&&faidaBetaMap_it->first<ts2){
            corrts = (Long64_t) faidaBetaMap_it->first;
            AIDASimpleStruct* beta = faidaBetaMap_it->second;
            double betax= beta->GetHitPositionX();
            double betay= beta->GetHitPositionY();
            short betaz= beta->GetHitPositionZ();
            if (corrts!=check_time&&betaz==impz){
                if (!((betax-impx>=-fmaxnpixels)&&(betax-impx<=fmaxnpixels)&&(betay-impy>=-fmaxnpixels)&&(betay-impy<=fmaxnpixels))){
                    faidaBetaMap_it++;
                    continue;
                }
                //! fill beta here
                //fh1->Fill((Long64_t)corrts-(Long64_t)ts);


                //flocalbeta->Clear();
                //flocalbeta->CopyFromAIDA(beta);
                //ftree->Fill();

                IonBeta* simplestructbeta=(IonBeta*)flocalbeta->ConstructedAt(nionbetacorr);
                simplestructbeta->CopyFromAIDA(beta);

                nionbetacorr++;
            }
            faidaBetaMap_it++;
        }
        if (nionbetacorr>0) ncorrwbeta++;
        ftree->Fill();
        k++;
    }
}

void Merger::DoMergeSingle()
{
    Int_t k=0;
    Int_t ktotal=faidaIonMap.size();
    Int_t ncorrwbigrips=0;
    Int_t ncorrwneutron=0;
    Int_t ncorrwbeta=0;


    Long64_t deadtime_start=-9999;
    Long64_t deadtime_reset=-9999;
    ff11vetodeadtime=0;
    ff11vetototaltime = 0;
    unsigned long long lastts = 0;
    for (fF11MapR_it=fF11RMap.begin();fF11MapR_it!=fF11RMap.end();fF11MapR_it++){
        unsigned long long ts=fF11MapR_it->first;
        unsigned long long entry=fF11MapR_it->second;
        if (ff11vetototaltime==0) ff11vetototaltime=(Long64_t)ts;
        lastts=ts;
        //ftrAnc->GetEntry(entry);
        //cout<<"ts="<<ts<<"\te="<<fanc->GetEnergy()<<endl;
        //! estimate dead time using paralyzable model(note: fNeuAncTWlow must be 0!)
        if (deadtime_reset<ts){
            if(deadtime_start>0) {
                ff11vetodeadtime+=deadtime_reset-deadtime_start;//prev window
                //cout<<"deadtime=\t"<<deadtime_reset-deadtime_start<<endl;
            }
            deadtime_start=ts;
            deadtime_reset=ts+fNeuAncTWup;
        }else{
            deadtime_reset=ts+fNeuAncTWup;
        }
        //unsigned int entry=fF11MapR_it->second;
        //ftrAnc->GetEntry(entry);
        //if (fanc->GetEnergy()<fmineneuvetodown) continue;
        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetF11Time();
                if (currtdiff<999998.) neuhit->SetF11Time(tdiff);
                else if (tdiff<currtdiff) neuhit->SetF11Time(tdiff);
                ncorr++;
            }
            fhe3Map_it++;
        }

        //! add all veto map
        fvetoMap.insert(make_pair(ts,entry));
    }
    ff11vetototaltime=(Long64_t)lastts+fNeuAncTWup-ff11vetototaltime;
    ff11vetodeadtime+=fNeuAncTWup;

    cout<<"F11R veto: deadtime="<<ff11vetodeadtime<<"\t---total time---\t"<<ff11vetototaltime<<endl;

    deadtime_start=-9999;
    deadtime_reset=-9999;
    fdownstreamvetodeadtime=0;
    fdownstreamvetototaltime=0;
    lastts=0;
    //! Neutron correlation with downstream veto
    for (fVetoDownMap_it=fVetoDownMap.begin();fVetoDownMap_it!=fVetoDownMap.end();fVetoDownMap_it++){
        unsigned long long ts=fVetoDownMap_it->first;
        unsigned int entry=fVetoDownMap_it->second;
        if (fdownstreamvetototaltime==0) fdownstreamvetototaltime=(Long64_t)ts;
        lastts=ts;
        ftrAnc->GetEntry(entry);
        if (fanc->GetEnergy()<fmineneuvetodown) continue;

        if (deadtime_reset<ts){
            if(deadtime_start>0) {
                fdownstreamvetodeadtime+=deadtime_reset-deadtime_start;//prev window
                //cout<<"deadtime=\t"<<deadtime_reset-deadtime_start<<endl;
            }
            deadtime_start=ts;
            deadtime_reset=ts+fNeuAncTWup;
        }else{
            deadtime_reset=ts+fNeuAncTWup;
        }
        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetDownstreamVetoTime();
                if (currtdiff<999998.) neuhit->SetDownstreamVetoTime(tdiff);
                else if (tdiff<currtdiff) neuhit->SetDownstreamVetoTime(tdiff);
                ncorr++;
                fh1->Fill(tdiff);
            }
            fhe3Map_it++;
        }
        //! add all veto map
        fvetoMap.insert(make_pair(ts,entry));
    }
    fdownstreamvetototaltime=(Long64_t)lastts+fNeuAncTWup-fdownstreamvetototaltime;
    fdownstreamvetodeadtime+=fNeuAncTWup;
    cout<<"Downstream veto:deadtime="<<fdownstreamvetodeadtime<<"\t---total time---\t"<<fdownstreamvetototaltime<<endl;

    //! allveto map
    deadtime_start=-9999;
    deadtime_reset=-9999;
    fvetodeadtime=0;
    fvetototaltime=0;
    lastts=0;
    //! Neutron correlation with downstream veto
    for (fvetoMap_it=fvetoMap.begin();fvetoMap_it!=fvetoMap.end();fvetoMap_it++){
        unsigned long long ts=fvetoMap_it->first;
        //unsigned int entry=fvetoMap_it->second;
        if (fvetototaltime==0) fvetototaltime=(Long64_t)ts;
        lastts=ts;
        //ftrAnc->GetEntry(entry);

        if (deadtime_reset<ts){
            if(deadtime_start>0) {
                fvetodeadtime+=deadtime_reset-deadtime_start;//prev window
                //cout<<"deadtime=\t"<<deadtime_reset-deadtime_start<<endl;
            }
            deadtime_start=ts;
            deadtime_reset=ts+fNeuAncTWup;
        }else{
            deadtime_reset=ts+fNeuAncTWup;
        }
        //! Correlate f11 with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetFinalVetoTime();
                if (currtdiff<999998.) neuhit->SetFinalVetoTime(tdiff);
                else if (tdiff<currtdiff) neuhit->SetFinalVetoTime(tdiff);
                ncorr++;
                fh1->Fill(tdiff);
            }
            fhe3Map_it++;
        }
    }
    fvetototaltime=(Long64_t)lastts+fNeuAncTWup-fvetototaltime;
    fvetodeadtime+=fNeuAncTWup;
    cout<<"Final veto: deadtime="<<fvetodeadtime<<"\t---total time---\t"<<fvetototaltime<<endl;

    //! Counter for neutron veto
    Int_t nvetoneu=0;
    Int_t ntotalneu=0;
    for (fhe3Map_it=fhe3Map.begin();fhe3Map_it!=fhe3Map.end();fhe3Map_it++){
        BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
        if (neuhit->GetEnergy()<860&&neuhit->GetEnergy()>160) {
            if (neuhit->GetFinalVetoTime()>0) nvetoneu++;
            ntotalneu++;
        }
    }
    cout<<"Total number of neutron="<<ntotalneu<<" ,vetoed neutrons="<<nvetoneu<<endl;


    //! AIDA-BIGRIPS correlation
    for (faidaIonMap_it=faidaIonMap.begin();faidaIonMap_it!=faidaIonMap.end();faidaIonMap_it++){
        flocalimp->Clear();
        flocalbeta->Clear("C");
        //flocalbeta->Clear();
        if (k%1000==0) cout<<k<<"/"<<ktotal<<"\t ncorr w bigrips "<<ncorrwbigrips<<endl;
        //if (k>500) break;
        unsigned long long ts=faidaIonMap_it->first;
        AIDASimpleStruct* ion=(AIDASimpleStruct*) faidaIonMap_it->second;

        ImplantCorrelationVector* corrvector=new ImplantCorrelationVector;
        corrvector->correntrybrips=-9999;
        corrvector->correntryf11r=-9999;
        corrvector->correntrydEbot=-9999;
        corrvector->correntrydEtop=-9999;
        corrvector->correntryvetodown=-9999;

        //! Correlate imp with bigrips
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        fbigripsMap_it = fbigripsMap.lower_bound(ts1);
        while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
            corrts = (Long64_t) fbigripsMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntrybrips = (int) fbigripsMap_it->second;

                //! fill data here
                //ftrBigrips->GetEvent(correntry);
                ncorr++;
                break;
            }
            fbigripsMap_it++;
        }

        //! Correlate imp with f11r
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryf11r = (int) fF11MapR_it->second;
                ncorr++;
                break;
            }
            fF11MapR_it++;
        }

        //! Correlate imp with topde
        ts1 = (Long64_t)ts - (Long64_t)fIondETWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIondETWup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fdETopMap_it = fdETopMap.lower_bound(ts1);
        while(fdETopMap_it!=fdETopMap.end()&&fdETopMap_it->first<ts2){
            corrts = (Long64_t) fdETopMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntrydEtop = (int) fdETopMap_it->second;
                ncorr++;
                break;
            }
            fdETopMap_it++;
        }

        //! Correlate imp with botde
        ts1 = (Long64_t)ts - (Long64_t)fIondETWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIondETWup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fdEBotMap_it = fdEBotMap.lower_bound(ts1);
        while(fdEBotMap_it!=fdEBotMap.end()&&fdEBotMap_it->first<ts2){
            corrts = (Long64_t) fdEBotMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntrydEbot = (int) fdEBotMap_it->second;
                ncorr++;
                break;
            }
            fdEBotMap_it++;
        }

        //! Correlate imp with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryvetodown = (int) fVetoDownMap_it->second;
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }

        faidaImplantMap.insert(make_pair(ts,make_pair(corrvector,ion)));
        if (ncorr>0) ncorrwbigrips++;
        k++;
    }


    //! Build Decay
    //!**************
    ktotal=faidaBetaMap.size();
    k=0;
    for (faidaBetaMap_it=faidaBetaMap.begin();faidaBetaMap_it!=faidaBetaMap.end();faidaBetaMap_it++){
        flocalbetaS->Clear();
        if (k%10000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorrwbeta<<"\t ncorr w neutron "<<ncorrwneutron<<endl;
        //if (k>500) break;
        unsigned long long ts=faidaBetaMap_it->first;
        AIDASimpleStruct* beta=(AIDASimpleStruct*) faidaBetaMap_it->second;
        double betax=beta->GetHitPositionX();
        double betay=beta->GetHitPositionY();
        short betaz=beta->GetHitPositionZ();//! correct dZ


        //! Time correlation with neutron
        int ndelayedneutron=0;
        //! with neutron
        Long64_t ts1 = ts - fNeuBetaTWup;
        Long64_t ts2 = ts + fNeuBetaTWlow;// this is upper :)
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                BELENHit* neuhitclone=new BELENHit;
                neuhit->Copy(*neuhitclone);
                if ((corrts-(Long64_t)ts)>=fNeuBetaoffset) flocalbetaS->AddNeutronForward(neuhitclone);
                else flocalbetaS->AddNeutronBackward(neuhitclone);
                //fh1->Fill(neutronts-(Long64_t)ts);
                ndelayedneutron++;
            }
            fhe3Map_it++;
        }
        if (ndelayedneutron>0) ncorrwneutron++;

        //!******************************ANCINARY detector****************
        //! Correlate imp with f11r
        ts1 = (Long64_t)ts - (Long64_t)fBetaAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaAncTWup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts = (Long64_t) fF11MapR_it->first;
            correntry = fF11MapR_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                ftrAnc->GetEvent(correntry);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalbetaS->AddAnc(hit);
                ncorr++;
                break;
            }
            fF11MapR_it++;
        }

        //! Correlate with vetodown
        ts1 = (Long64_t)ts - (Long64_t)fBetaAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaAncTWup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts = (Long64_t) fVetoDownMap_it->first;
            correntry = fVetoDownMap_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                ftrAnc->GetEvent(correntry);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalbetaS->AddAnc(hit);
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }


        //!******************************Implantation****************
        //! Time correlation with implantation (Build Decay Curve)
        ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;

        //! Correlate with beta
        ts1 = (Long64_t)ts - fIonBetaTWlow;
        ts2 = (Long64_t)ts + fIonBetaTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        nionbetacorr=0;
        faidaImplantMap_it = faidaImplantMap.lower_bound(ts1);
        while(faidaImplantMap_it!=faidaImplantMap.end()&&faidaImplantMap_it->first<ts2){
            flocalimp->Clear();
            corrts = (Long64_t) faidaImplantMap_it->first;
            AIDASimpleStruct* imp = faidaImplantMap_it->second.second;
            ImplantCorrelationVector* corrvectorimp= (ImplantCorrelationVector*)faidaImplantMap_it->second.first;

            double impx= imp->GetHitPositionX();
            double impy= imp->GetHitPositionY();
            short impz= imp->GetHitPositionZ();
            if (corrts!=check_time&&betaz==impz){
                if (!((betax-impx>=-fmaxnpixels)&&(betax-impx<=fmaxnpixels)&&(betay-impy>=-fmaxnpixels)&&(betay-impy<=fmaxnpixels))){
                    faidaImplantMap_it++;
                    continue;
                }
                check_time=corrts;
                //! fill bigrips data here
                if (corrvectorimp->correntrybrips>=0) {
                    ftrBigrips->GetEvent(corrvectorimp->correntrybrips);
                    flocalimp->AddBeam(*fbigrips);
                }

                //! fill anc data here
                if (corrvectorimp->correntryf11r>=0){
                    ftrAnc->GetEvent(corrvectorimp->correntryf11r);
                    BELENHit* hit=new BELENHit;
                    fanc->Copy(*hit);
                    flocalimp->AddAnc(hit);
                }
                if (corrvectorimp->correntrydEtop>=0){
                    ftrAnc->GetEvent(corrvectorimp->correntrydEtop);
                    BELENHit* hit=new BELENHit;
                    fanc->Copy(*hit);
                    flocalimp->AddAnc(hit);
                }
                if (corrvectorimp->correntrydEbot>=0){
                    ftrAnc->GetEvent(corrvectorimp->correntrydEbot);
                    BELENHit* hit=new BELENHit;
                    fanc->Copy(*hit);
                    flocalimp->AddAnc(hit);
                }
                if (corrvectorimp->correntryvetodown>=0){
                    ftrAnc->GetEvent(corrvectorimp->correntryvetodown);
                    BELENHit* hit=new BELENHit;
                    fanc->Copy(*hit);
                    flocalimp->AddAnc(hit);
                }


                //! fill beta data here

                flocalbetaS->CopyFromAIDA(beta);
                flocalimp->CopyFromAIDA(imp);               
                ftree->Fill();

                //fh1->Fill((Long64_t)ts-(Long64_t)corrts);

                nionbetacorr++;
            }
            faidaImplantMap_it++;
        }
        if (nionbetacorr>0) ncorrwbeta++;
        k++;
    }

}


void Merger::DoMergeTest()
{
    Int_t k=0;
    Int_t ktotal=faidaBetaMap.size();
    Int_t ncorrwbigrips=0;
    Int_t ncorrwbeta=0;
    for (faidaBetaMap_it=faidaBetaMap.begin();faidaBetaMap_it!=faidaBetaMap.end();faidaBetaMap_it++){
        flocalimp->Clear();
        flocalbeta->Clear("C");
        //flocalbeta->Clear();
        if (k%1000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorrwbeta<<"\t ncorr w bigrips "<<ncorrwbigrips<<endl;
        //if (k>500) break;
        unsigned long long ts=faidaBetaMap_it->first;
        AIDASimpleStruct* ion=(AIDASimpleStruct*) faidaBetaMap_it->second;
        double impx=ion->GetHitPositionX();
        double impy=ion->GetHitPositionY();
        short impz=ion->GetHitPositionZ()+ion->GetDZ();//! correct dZ

        Long64_t ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;

        //! Correlate with beta
        ts1 = (Long64_t)ts - fIonBetaTWlow;
        ts2 = (Long64_t)ts + fIonBetaTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        nionbetacorr=0;
        faidaIonMap_it = faidaIonMap.lower_bound(ts1);
        while(faidaIonMap_it!=faidaIonMap.end()&&faidaIonMap_it->first<ts2){
            corrts = (Long64_t) faidaIonMap_it->first;
            AIDASimpleStruct* beta = faidaIonMap_it->second;
            double betax= beta->GetHitPositionX();
            double betay= beta->GetHitPositionY();
            short betaz= beta->GetHitPositionZ();
            if (corrts!=check_time&&betaz==impz){
                if (!((betax-impx>=-fmaxnpixels)&&(betax-impx<=fmaxnpixels)&&(betay-impy>=-fmaxnpixels)&&(betay-impy<=fmaxnpixels))){
                    faidaIonMap_it++;
                    continue;
                }
                check_time=corrts;
                //! fill beta here
                //fh1->Fill((Long64_t)ts-(Long64_t)corrts);
                nionbetacorr++;
            }
            faidaIonMap_it++;
        }
        if (nionbetacorr>0) ncorrwbeta++;
        k++;
    }
}


void Merger::DoMergeYOnly()
{
    Int_t k=0;
    Int_t ktotal=faidaIonMap.size();
    Int_t ncorrwbigrips=0;
    Int_t ncorrwbeta=0;
    for (faidaIonMap_it=faidaIonMap.begin();faidaIonMap_it!=faidaIonMap.end();faidaIonMap_it++){
        flocalimp->Clear();
        flocalbeta->Clear("C");
        //flocalbeta->Clear();
        if (k%1000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorrwbeta<<"\t ncorr w bigrips "<<ncorrwbigrips<<endl;
        //if (k>500) break;
        unsigned long long ts=faidaIonMap_it->first;
        AIDASimpleStruct* ion=(AIDASimpleStruct*) faidaIonMap_it->second;
        double impx=ion->GetHitPositionX();
        double impy=ion->GetHitPositionY();
        short impz=ion->GetHitPositionZ()+ion->GetDZ();//! correct dZ

        //! fill imp data here
        //AIDASimpleStruct* iontofill = flocalimp->GetIonBetaCluster();
        //ion->Copy(*iontofill);
        flocalimp->CopyFromAIDA(ion);
        //! Correlate imp with bigrips
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fIonPidTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fIonPidTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fbigripsMap_it = fbigripsMap.lower_bound(ts1);
        while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
            corrts = (Long64_t) fbigripsMap_it->first;
            correntry = fbigripsMap_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                //! fill data here
                ftrBigrips->GetEvent(correntry);
                flocalimp->AddBeam(*fbigrips);
                ncorr++;
                break;
            }
            fbigripsMap_it++;
        }
        if (ncorr>0) ncorrwbigrips++;

        //! Correlate with beta
        ts1 = (Long64_t)ts - fIonBetaTWlow;
        ts2 = (Long64_t)ts + fIonBetaTWup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        nionbetacorr=0;
        faidaBetaMap_it = faidaBetaMap.lower_bound(ts1);
        while(faidaBetaMap_it!=faidaBetaMap.end()&&faidaBetaMap_it->first<ts2){
            corrts = (Long64_t) faidaBetaMap_it->first;
            AIDASimpleStruct* beta = faidaBetaMap_it->second;
            double betax= beta->GetHitPositionX();
            double betay= beta->GetHitPositionY();
            short betaz= beta->GetHitPositionZ();
            if (corrts!=check_time&&betaz==impz){
                if (!((betay-impy>=-fmaxnpixels)&&(betay-impy<=fmaxnpixels))){
                    faidaBetaMap_it++;
                    continue;
                }
                //! fill beta here
                fh1->Fill((Long64_t)corrts-(Long64_t)ts);

                //flocalbeta->Clear();
                //flocalbeta->CopyFromAIDA(beta);
                //ftree->Fill();

                IonBeta* simplestructbeta=(IonBeta*)flocalbeta->ConstructedAt(nionbetacorr);
                simplestructbeta->CopyFromAIDA(beta);

                nionbetacorr++;
            }
            faidaBetaMap_it++;
        }
        if (nionbetacorr>0) ncorrwbeta++;
        ftree->Fill();
        k++;
    }
}

void Merger::DoSeparatePID()
{
    for (Long64_t jentry=0;jentry<fnentriesMerged;jentry++){        
        ftrMerged->GetEntry(jentry);
        double zet=flocalimp->GetBeamHit()->zet;
        double aoq=flocalimp->GetBeamHit()->aoq;
        for (Int_t j=0;j<nri;j++){
            if (!enablepid2[j]) continue;
            if (cutg[j]->IsInside(aoq,zet)){
                ftreeRI[j]->Fill();
            }
        }
        ftree->Fill();
    }
}
