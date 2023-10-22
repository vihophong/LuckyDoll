#include "Merger.h"
Merger::Merger():fbigrips(),decay()
{
    fmaxmult=400;
    fmaxnpixels=4.;
    //fIonBetaTWlow=10000000000;
    fIonBetaTWlow=20000000000;
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

    fNeuAncTWup = 200000; //200 us
    fNeuAncTWlow = 0;

    //fNeuAncTWup = 1000000;
    //fNeuAncTWlow = 1000000;

    fNeuBetaTWup = 200000; //200 us
    fNeuBetaTWlow = 200000; //200 us
    fNeuBetaoffset = -22000;
    //fNeuBetaoffset = 0;//for simulation

    fBetaGammaTWup = 10000;
    fBetaGammaTWlow = 40000;
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


    //! final cut
    dtioncut=600;
    sumexyrankcut=0;
    xytdiffcut=5000;
    lightp_nzcut=6;
    neutronecut[0]=160;
    neutronecut[1]=860;

    deltaxcut=5.;
    deltaycut=5.;


    //! gamma calibration
    Double_t sep[8] = {4.912531E+05,4.998953E+05,5.248885E+05,7.622378E+05,1.108068E+06,1.147523E+06,1.164587E+06,9.553438E+05}; //separation points
    Double_t low_offset[8] = {0.310056,0.602804	,0.594689,0.00575245,-0.409776,-0.739931,-0.231043,-0.410475};//low offset
    Double_t low_gain[8] = {0.00158189,0.00154124,0.00153239,0.000668465,0.000707824,0.000687903,0.000672509,0.000700449};//low gain
    Double_t low_se[8]= {2.5638600000E-12,5.1029500000E-12,8.0936100000E-12,-2.5440800000E-12,-3.2340000000E-12,-4.4067300000E-12,-3.9413800000E-12,-2.6372600000E-12};// low second order
    Double_t high_offset[8] = {1.32609,4.65127,-17.4998,4.21305,-1.68140e+00,-1.13557,-2.79797,6.44019};//high offset
    Double_t high_gain[8] = {0.00158113,0.00153239,0.00159751,0.000659638,7.05388e-04,0.000682542,0.000670503,0.000688148};//high gain
    Double_t high_se[8] = {-9.9227100000E-14,6.6060100000E-12,-5.0294000000E-11,1.7949000000E-12,1.00018e-16,5.6552400000E-13,-3.2623400000E-13,2.7326500000E-12};//high second order
    Double_t gain_old[8] =
    {0.0015836325,0.0015439776,0.0015329461,0.000664755,0.00070,0.000683305,0.0006685615,0.0006968099};
    Double_t offset_old[8] =
    {-0.6017825233,-0.2629124127,0.2931096096,-0.0427639699,-0.5577366054,-0.6254548603,-0.1966707051,-0.2212697416};


    for (int i=0;i<8;i++){
        fsep[i]=sep[i];
        flow_offset[i]=low_offset[i];
        flow_gain[i]=low_gain[i];
        flow_se[i]=low_se[i];
        fhigh_offset[i]=high_offset[i];
        fhigh_gain[i]=high_gain[i];
        fhigh_se[i]=high_se[i];
        fcgainold[i]=gain_old[i];
        fcoffsetold[i]=offset_old[i];
    }

    //! stuff for rejecting noise events associated with implantation
    fimpnoisefilter_dxy = 4;//3 pixel from implantation/punching through position
    fimpnoisefilter_dz = 0;//1 dssd downstream of implantation position
    fimpnoisefilter_dt = 50000000; //unit ns: reject 50 usecond
}

Merger::~Merger()
{
    delete fh1;
    delete fh2;
    delete fh2d1;
    delete fh2d2;
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
        ftreeimplantRI[i]=new TTree(Form("treeimp%s",(char*)nameri[i].Data()),Form("treeimp%s",(char*)nameri[i].Data()));
    }
}

void Merger::InitPIDSep()
{
    ResetSimpleData();

    flocalimparray = new TClonesArray("IonBetaMult");
    flocalimparray->Clear("C");
    flocalbetaS = new IonBeta;
    flocalbetaS->Clear();
    //! init aida
    fmergedFile = new TFile(finputMerged);
    fmergedFile->GetObject("tree",ftrMerged);
    ftrMerged->SetBranchAddress("idx",&nionbetacorr);
    ftrMerged->SetBranchAddress("ion",&flocalimparray);
    ftrMerged->GetBranch("beta")->SetAutoDelete(kFALSE);
    ftrMerged->SetBranchAddress("beta",&flocalbetaS);

    fnentriesMerged = ftrMerged->GetEntries();
    cout<<"Reading "<<fnentriesMerged<<" enetries in Merged tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftrMerged->GetEvent(i);
        Int_t nimpevt=flocalimparray->GetEntriesFast();
        for (Int_t j=0;j<nimpevt;j++){
            IonBetaMult* imp=(IonBetaMult*)flocalimparray->At(j);
            cout<<"ionts="<<imp->GetTimestamp()<<endl;
            cout<<"betats= "<<flocalbetaS->GetTimestamp()<<endl;
        }
    }

    //! init tree implant
    fmergedFile->GetObject("treeimplant",ftreeImplant);
    ftreeImplant->SetBranchAddress("implant",&flocalimp);
    fnentriesImp=ftreeImplant->GetEntries();
    cout<<"Reading "<<fnentriesImp<<" enetries in Implant tree"<<endl;
    cout<<"Printing first few timestamp:"<<endl;
    for (Long64_t i=0;i<10;i++){
        ftreeImplant->GetEvent(i);
        cout<<"ionts= "<<flocalimp->GetTimestamp()<<endl;
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
    fh1=new TH1F("h1","h1",2000,-10e9,10e9);
    fh2=new TH1F("h2","h2",2000,-10e9,10e9);
}

void Merger::Init()
{
    ftree = 0;
    ftreeImplant = 0;
    ftreeBeta = 0;
    ftreeNeutron = 0;



    flocalimparray = new TClonesArray("IonBetaMult");
    flocalimparray->Clear("C");

    flocalbetaS = new IonBeta;
    flocalbetaS->Clear();


    flocalimp = new IonBetaMult;
    flocalimp->Clear();


    flocalneutron = new BELENHit;

    nionbetacorr=0;

    faida = NULL;
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


    //! stuffs for deadtime estimation using pulser

    ftsbeginpulser=0;
    ftsendpulser=0;
    for (Int_t i=0;i<141;i++){
        fhpulser[i]=new TH1F(Form("hpulser%d",i+1),Form("hpulser%d",i+1),4000,1500,3500);
        fhpulserall[i]=new TH1F(Form("fhpulserall%d",i+1),Form("fhpulserall%d",i+1),4000,1500,3500);
    }
    fh2deadtime=new TH2F("h2deadtime","h2deadtime with veto using calculated start-stop time",141,0,141,2000,0,20);
    fh1deadtime=new TH1F("h1deadtime","h1deadtime with veto using calculated start-stop time",2000,0,20);

    fh2deadtime2=new TH2F("h2deadtime2","h2deadtime without veto using calculated start-stop time",141,0,141,2000,0,20);
    fh1deadtime2=new TH1F("h1deadtime2","h1deadtime without veto using calculated start-stop time",2000,0,20);

    fh2deadtime3=new TH2F("h2deadtime3","h2deadtime with veto using dtpulser",141,0,141,2000,0,20);
    fh1deadtime3=new TH1F("h1deadtime3","h1deadtime with veto using dtpulser",2000,0,20);

    fh2deadtime4=new TH2F("h2deadtime4","h2deadtime without veto using dtpulser",141,0,141,2000,0,20);
    fh1deadtime4=new TH1F("h1deadtime4","h1deadtime without veto using dtpulser",2000,0,20);


    fh1dtpulser=new TH1F("h1dtpulser","h1dtpulser",2000,0,2000);
    fh1=new TH1F("h1","h1",2000,0,1e9);
    fh2=new TH1F("h2","h2",1000,0,600);
    fh2d1=new TH2F("h2d1","h2d1",1000,0,100,1000,0,100);
    fh2d2=new TH2F("h2d2","h2d2",1000,0,400,1000,0,100);

    fNeuBetaTWup= fNeuBetaTWup-fNeuBetaoffset;
    fNeuBetaTWlow= fNeuBetaTWlow+fNeuBetaoffset;//this is upper :)
    cout<<fNeuBetaTWup<<"\t"<<fNeuBetaoffset<<"\t"<<fNeuBetaTWlow<<endl;

}



void Merger::BookTreeSingle(TTree* tree){
    ftree = tree;
    ftree->Branch("idx",&nionbetacorr,"idx/I");
    ftree->Branch("ion",&flocalimparray);
    ftree->Branch("beta",&flocalbetaS);
    ftree->BranchRef();
}
void Merger::BookTreeNeutron(TTree* tree)
{
    ftreeNeutron = tree;
    ftreeNeutron->Branch("neutron",&flocalneutron);
    ftreeNeutron->BranchRef();
}
void Merger::BookTreeImplant(TTree* tree)
{
    ftreeImplant = tree;
    ftreeImplant->Branch("implant",&flocalimp);
    ftreeImplant->BranchRef();
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

void Merger::BookPIDSepSimpleTree(){
    for (Int_t i=0;i<nri;i++){
        ftreeRI[i]->Branch("idx",&nionbetacorr,"idx/I");
        ftreeRI[i]->Branch("ion",&flocalimp);
        ftreeRI[i]->Branch("beta",&flocalbetaS);
        ftreeRI[i]->Branch("decay",&decay,"evt/l:ts/l:t/D:x/D:y/D:ex/D:ey/D:ion_x/D:ion_y/D:ion_ex/D:ion_ey/D:zet/D:aoq/D:deltaxy/D:z/S:ion_z/S:multx/S:multy/S:multz/S:ndecay/S:nbeta/S");
        ftreeRI[i]->Branch("gc_hit",&decay.gc_hit,"gc_hit/I");
        ftreeRI[i]->Branch("gc_E",decay.gc_E,"gc_E[gc_hit]/D");
        ftreeRI[i]->Branch("gc_T",decay.gc_T,"gc_T[gc_hit]/D");
        ftreeRI[i]->Branch("gc_ch",decay.gc_ch,"gc_ch[gc_hit]/I");
        ftreeRI[i]->Branch("neu_hit",&decay.neu_hit,"neu_hit/I");
        ftreeRI[i]->Branch("neu_E",decay.neu_E,"neu_E[neu_hit]/D");
        ftreeRI[i]->Branch("neu_T",decay.neu_T,"neu_T[neu_hit]/D");
        ftreeRI[i]->Branch("neu_ch",decay.neu_ch,"neu_ch[neu_hit]/I");
        ftreeRI[i]->Branch("neu_x",decay.neu_x,"neu_x[neu_hit]/D");
        ftreeRI[i]->Branch("neu_y",decay.neu_y,"neu_y[neu_hit]/D");
        ftreeRI[i]->Branch("neub_hit",&decay.neub_hit,"neub_hit/I");
        //ftreeRI[i]->BranchRef();
        ftreeimplantRI[i]->Branch("implant",&flocalimp);
    }
    ftreeimplantAll=new TTree("treeimpall","treeimpall");
    ftreeimplantAll->Branch("implant",&flocalimp);

    ftree=new TTree("tree","tree");
    ftree->Branch("decay",&decay,"evt/l:ts/l:t/D:x/D:y/D:ex/D:ey/D:ion_x/D:ion_y/D:ion_ex/D:ion_ey/D:zet/D:aoq/D:deltaxy/D:z/S:ion_z/S:multx/S:multy/S:multz/S:ndecay/S:nbeta/S");
    ftree->Branch("gc_hit",&decay.gc_hit,"gc_hit/I");
    ftree->Branch("gc_E",decay.gc_E,"gc_E[gc_hit]/D");
    ftree->Branch("gc_T",decay.gc_T,"gc_T[gc_hit]/D");
    ftree->Branch("gc_ch",decay.gc_ch,"gc_ch[gc_hit]/I");
    ftree->Branch("neu_hit",&decay.neu_hit,"neu_hit/I");
    ftree->Branch("neu_E",decay.neu_E,"neu_E[neu_hit]/D");
    ftree->Branch("neu_T",decay.neu_T,"neu_T[neu_hit]/D");
    ftree->Branch("neu_ch",decay.neu_ch,"neu_ch[neu_hit]/I");
    ftree->Branch("neu_x",decay.neu_x,"neu_x[neu_hit]/D");
    ftree->Branch("neu_y",decay.neu_y,"neu_y[neu_hit]/D");
    ftree->Branch("neub_hit",&decay.neub_hit,"neub_hit/I");
}


void Merger::ResetSimpleData(){
    decay.evt=0;
    decay.ts=0;
    decay.t=-9999.;decay.x=-9999.;decay.y=-9999.;decay.ex=-9999.;decay.ey=-9999.;decay.ion_x=-9999.;decay.ion_y=-9999.;decay.ion_ex=-9999.;decay.ion_ey=-9999.;
    decay.zet=-9999.;decay.aoq=-9999.;decay.deltaxy=-9999.;
    decay.z=-9999;decay.ion_z=-9999;decay.multx=0;decay.multy=0;decay.multz=0;decay.ndecay=-9999;decay.nbeta=-9999;

    decay.gc_hit=0;

    for (int i=0;i<kMaxGamma;i++){
        decay.gc_E[i]=-9999.;
        decay.gc_T[i]=-9999.;
        decay.gc_ch[i]=-9999;
    }

    decay.neu_hit=0;
    for (int i=0;i<kMaxNeutron;i++){
        decay.neu_E[i]=-9999.;
        decay.neu_T[i]=-9999.;
        decay.neu_ch[i]=-9999;
        decay.neu_x[i]=-9999.;
        decay.neu_y[i]=-9999.;
    }
    decay.neub_hit=0;

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

    //! stuff for calculating cycle of time stamp reset
    unsigned long long prevts=0;
    unsigned short currcycle=0;
    for (unsigned int jentry = sstartN;jentry < sstopN;jentry++){
        ftrNeutron->GetEvent(jentry);
        if (fneutron->GetTimestamp()<prevts) currcycle++;// if time jump detected -> increase cycle by 1
        prevts=fneutron->GetTimestamp();
        BELENHit* hit=new BELENHit;        
        fneutron->Copy(*hit);
        hit->SetHitsAdded(currcycle);
        fhe3Map.insert(make_pair(fneutron->GetTimestamp(),hit));
        //if (fneutron->GetEnergy()<fmaxneue&&fneutron->GetEnergy()>fminneue) fhe3Map.insert(make_pair(fneutron->GetTimestamp(),hit));
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

    prevts=0;
    currcycle=0;    
    Long64_t prevtsdtpulser=0;
    Int_t ncountsDTpulser = 0;
    for (unsigned int jentry = sstartA;jentry < sstopA;jentry++){
        ftrAnc->GetEvent(jentry);

        if (fanc->GetTimestamp()<prevts) currcycle++;// if time jump detected -> increase cycle by 1
        prevts=fanc->GetTimestamp();
        //fancMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==1) {
            fF11RMap.insert(make_pair(fanc->GetTimestamp(),jentry));            
        }
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==2) fF11LMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==1) fF11LRMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==2&&fanc->GetID()==1) fVetoTopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==2) fVetoBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==3&&fanc->GetID()==1) fdETopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==3&&fanc->GetID()==2) fdEBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==4) fVetoDownMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==5) {
            if (ftsbeginpulser==0&&ncountsDTpulser==20) ftsbeginpulser=fanc->GetTimestamp();
            ftsendpulser=fanc->GetTimestamp();

            BELENHit* hit=new BELENHit;
            fanc->Copy(*hit);
            hit->SetID(141);
            hit->SetHitsAdded(currcycle);
            fdtpulserMap.insert(make_pair(fanc->GetTimestamp(),hit));

            if (ftsbeginpulser>0 && prevtsdtpulser!=0) {
                fh1dtpulser->Fill(fanc->GetEnergy());//for dead time calculation
                fh1->Fill(ftsendpulser-prevtsdtpulser);
            }
            prevtsdtpulser=ftsendpulser;
            ncountsDTpulser++;
        }
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

    for (fF11MapLR_it=fF11LRMap.begin();fF11MapLR_it!=fF11LRMap.end();fF11MapLR_it++){
        unsigned long long ts=fF11MapLR_it->first;
        unsigned long long entry=fF11MapLR_it->second;
        if (ff11vetototaltime==0) ff11vetototaltime=(Long64_t)ts;
        lastts=ts;
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

        ftrAnc->GetEntry(entry);

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
            //if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetF11Time();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetF11Time(tdiff);
                //if (ncorr==0) fh2->Fill(fanc->GetEnergy());
                fh2->Fill(tdiff);

            //}
            fhe3Map_it++;
        }
        //fh1->Fill(ncorr);

        //! correlation of with deadtime pulser map
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fdtpulserMap_it = fdtpulserMap.lower_bound(ts1);
        while(fdtpulserMap_it!=fdtpulserMap.end()&&fdtpulserMap_it->first<ts2){
            corrts = (Long64_t) fdtpulserMap_it->first;
            BELENHit* neuhit = (BELENHit*) fdtpulserMap_it->second;
            //if (corrts!=check_time){//reject duplicate time stamp, disable for pulser events
                check_time=corrts;
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetFinalVetoTime();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetF11Time(tdiff);
                ncorr++;
            //}
            fdtpulserMap_it++;
        }
        //! add all veto map
        fvetoMap.insert(make_pair(ts,entry));
    }
    ff11vetototaltime=(Long64_t)lastts+fNeuAncTWup-ff11vetototaltime;
    ff11vetodeadtime+=fNeuAncTWup;

    cout<<"F11R veto: deadtime="<<ff11vetodeadtime<<"\t---total time---\t"<<ff11vetototaltime<<endl;
    cout<<"Found "<<fF11LRMap.size()<<"L-R or of F11 over "<<fF11RMap.size()<<" events"<<endl;

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
        if (deadtime_reset<ts){
            if(deadtime_start>0) {
                fdownstreamvetodeadtime+=deadtime_reset-deadtime_start;//prev window
            }
            deadtime_start=ts;
            deadtime_reset=ts+fNeuAncTWup;
        }else{
            deadtime_reset=ts+fNeuAncTWup;
        }
        //! Correlate  with neutron
        Long64_t ts1 = (Long64_t)ts - (Long64_t)fNeuAncTWlow;
        Long64_t ts2 = (Long64_t)ts + (Long64_t)fNeuAncTWup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        Long64_t check_time = 0;
        fhe3Map_it = fhe3Map.lower_bound(ts1);
        while(fhe3Map_it!=fhe3Map.end()&&fhe3Map_it->first<ts2){
            corrts = (Long64_t) fhe3Map_it->first;
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            //if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetDownstreamVetoTime();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetDownstreamVetoTime(tdiff);
                if (ncorr==0&&neuhit->GetEnergy()>neutronecut[0]&&neuhit->GetEnergy()<neutronecut[1]) fh2d2->Fill(tdiff,fanc->GetEnergy());
                if (neuhit->GetEnergy()>neutronecut[0]&&neuhit->GetEnergy()<neutronecut[1]) ncorr++;
                ncorr++;

            //}
            fhe3Map_it++;
        }

        //! correlation of with deadtime pulser map
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fdtpulserMap_it = fdtpulserMap.lower_bound(ts1);
        while(fdtpulserMap_it!=fdtpulserMap.end()&&fdtpulserMap_it->first<ts2){
            corrts = (Long64_t) fdtpulserMap_it->first;
            BELENHit* neuhit = (BELENHit*) fdtpulserMap_it->second;
            //if (corrts!=check_time){//reject duplicate time stamp, disable for pulser events
                check_time=corrts;
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetFinalVetoTime();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetDownstreamVetoTime(tdiff);
                ncorr++;
            //}
            fdtpulserMap_it++;
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
    //! Neutron correlation with all veto
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
            //if (corrts!=check_time){//reject duplicate time stamp, disable for pulser events
                check_time=corrts;
                //! fill neutron here
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetFinalVetoTime();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetFinalVetoTime(tdiff);
                ncorr++;
            //}
            fhe3Map_it++;
        }


        //! correlation of with deadtime pulser map
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fdtpulserMap_it = fdtpulserMap.lower_bound(ts1);
        while(fdtpulserMap_it!=fdtpulserMap.end()&&fdtpulserMap_it->first<ts2){
            corrts = (Long64_t) fdtpulserMap_it->first;
            BELENHit* neuhit = (BELENHit*) fdtpulserMap_it->second;
            //if (corrts!=check_time){//reject duplicate time stamp, disable for pulser events
                check_time=corrts;
                double tdiff=(Double_t)(corrts-(Long64_t)ts)/1000.;//in us
                double currtdiff=neuhit->GetFinalVetoTime();
                if (currtdiff<999998.||tdiff<currtdiff) neuhit->SetFinalVetoTime(tdiff);
                ncorr++;
            //}
            fdtpulserMap_it++;
        }

    }
    fvetototaltime=(Long64_t)lastts+fNeuAncTWup-fvetototaltime;
    fvetodeadtime+=fNeuAncTWup;
    cout<<"Final veto: deadtime="<<fvetodeadtime<<"\t---total time---\t"<<fvetototaltime<<endl;


    //! Counter for neutron veto and fill pulser histogram
    Int_t nvetoneu=0;
    Int_t ntotalneu=0;
    for (fhe3Map_it=fhe3Map.begin();fhe3Map_it!=fhe3Map.end();fhe3Map_it++){
        BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
        Int_t id=neuhit->GetID()-1;
        if (neuhit->GetTimestamp()>ftsbeginpulser&&neuhit->GetTimestamp()<ftsendpulser) fhpulserall[id]->Fill(neuhit->GetEnergy());
        //if (neuhit->GetF11Time()<0&&neuhit->GetTimestamp()>ftsbeginpulser&&neuhit->GetTimestamp()<ftsendpulser) fhpulser[id]->Fill(neuhit->GetEnergy());
        if (neuhit->GetFinalVetoTime()<0&&neuhit->GetTimestamp()>ftsbeginpulser&&neuhit->GetTimestamp()<ftsendpulser) fhpulser[id]->Fill(neuhit->GetEnergy());
        if (neuhit->GetEnergy()<neutronecut[1]&&neuhit->GetEnergy()>neutronecut[0]){
            //if (neuhit->GetF11Time()>=0) nvetoneu++;
            if (neuhit->GetFinalVetoTime()>=0) nvetoneu++;
            ntotalneu++;
        }
    }

    for (fdtpulserMap_it=fdtpulserMap.begin();fdtpulserMap_it!=fdtpulserMap.end();fdtpulserMap_it++){
        BELENHit* neuhit = (BELENHit*) fdtpulserMap_it->second;   
        Int_t id=neuhit->GetID()-1;
        if (id!=140) cout<<"error! "<<id<<endl;
        if (neuhit->GetTimestamp()>ftsbeginpulser&&neuhit->GetTimestamp()<ftsendpulser) fhpulserall[id]->Fill(neuhit->GetEnergy());
        if (neuhit->GetF11Time()<0&&neuhit->GetTimestamp()>ftsbeginpulser&&neuhit->GetTimestamp()<ftsendpulser) fhpulser[id]->Fill(neuhit->GetEnergy());
    }

    cout<<"************Total number of neutron="<<ntotalneu<<" ,vetoed neutrons="<<nvetoneu<<endl;

    //! Write out neutron tree
    if (ftreeNeutron!=0){
        cout<<"MAKING NEUTRON TREE"<<endl;
        for (fhe3Map_it=fhe3Map.begin();fhe3Map_it!=fhe3Map.end();fhe3Map_it++){
            flocalneutron->Clear();
            BELENHit* neuhit = (BELENHit*) fhe3Map_it->second;
            neuhit->Copy(*flocalneutron);
            ftreeNeutron->Fill();
        }        
        cout<<"add pulser data at the end"<<endl;
        for (fdtpulserMap_it=fdtpulserMap.begin();fdtpulserMap_it!=fdtpulserMap.end();fdtpulserMap_it++){
            flocalneutron->Clear();
            BELENHit* neuhit = (BELENHit*) fdtpulserMap_it->second;
            neuhit->Copy(*flocalneutron);
            ftreeNeutron->Fill();
        }
    }


    //! AIDA-BIGRIPS correlation
    for (faidaIonMap_it=faidaIonMap.begin();faidaIonMap_it!=faidaIonMap.end();faidaIonMap_it++){
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

        //! Correlate imp with f11lr
        ts1 = (Long64_t)ts - (Long64_t)fIonAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fIonAncTWup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fF11MapLR_it = fF11LRMap.lower_bound(ts1);
        while(fF11MapLR_it!=fF11LRMap.end()&&fF11MapLR_it->first<ts2){
            corrts = (Long64_t) fF11MapLR_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryf11r = (int) fF11MapLR_it->second;
                ncorr++;
                break;
            }
            fF11MapLR_it++;
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

    //! Write out implantation tree
    if (ftreeImplant!=0){
        cout<<"MAKING IMPLANT TREE"<<endl;
        for (faidaImplantMap_it=faidaImplantMap.begin();faidaImplantMap_it!=faidaImplantMap.end();faidaImplantMap_it++){
            flocalimp->Clear();
            unsigned long long tsimp=faidaImplantMap_it->first;
            AIDASimpleStruct* imp = faidaImplantMap_it->second.second;
            ImplantCorrelationVector* corrvectorimp= (ImplantCorrelationVector*)faidaImplantMap_it->second.first;

            //! make a copy of implant object and dum into vector
            IonBetaMult* impcpy=new IonBetaMult;

            //! fill bigrips data here
            if (corrvectorimp->correntrybrips>=0) {
                ftrBigrips->GetEvent(corrvectorimp->correntrybrips);
                flocalimp->AddBeam(*fbigrips);
                impcpy->AddBeam(*fbigrips);
            }
            //! fill anc data here
            if (corrvectorimp->correntryf11r>=0){
                ftrAnc->GetEvent(corrvectorimp->correntryf11r);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalimp->AddAnc(hit);
                BELENHit* hit2=new BELENHit;
                fanc->Copy(*hit2);
                impcpy->AddAnc(hit2);
            }
            if (corrvectorimp->correntrydEtop>=0){
                ftrAnc->GetEvent(corrvectorimp->correntrydEtop);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalimp->AddAnc(hit);
                BELENHit* hit2=new BELENHit;
                fanc->Copy(*hit2);
                impcpy->AddAnc(hit2);

            }
            if (corrvectorimp->correntrydEbot>=0){
                ftrAnc->GetEvent(corrvectorimp->correntrydEbot);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalimp->AddAnc(hit);
                BELENHit* hit2=new BELENHit;
                fanc->Copy(*hit2);
                impcpy->AddAnc(hit2);
            }
            if (corrvectorimp->correntryvetodown>=0){
                ftrAnc->GetEvent(corrvectorimp->correntryvetodown);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalimp->AddAnc(hit);
                BELENHit* hit2=new BELENHit;
                fanc->Copy(*hit2);
                impcpy->AddAnc(hit2);
            }
            flocalimp->CopyFromAIDA(imp);
            impcpy->CopyFromAIDA(imp);
            faidaImplantMapFull.insert(make_pair(tsimp,impcpy));
            ftreeImplant->Fill();
        }
    }

    //! Build Decay
    //!**************
    ktotal=faidaBetaMap.size();
    k=0;
    for (faidaBetaMap_it=faidaBetaMap.begin();faidaBetaMap_it!=faidaBetaMap.end();faidaBetaMap_it++){
        //! vector for ion-beta correlation
        flocalbetaS->Clear();
        flocalimparray->Clear("C");
        if (k%10000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorrwbeta<<"\t ncorr w neutron "<<ncorrwneutron<<endl;
        //if (k>1) break;
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
            if (!(neuhit->GetEnergy()<fmaxneue&&neuhit->GetEnergy()>fminneue))
            {
                fhe3Map_it++;
                continue;
            }
            if (corrts!=check_time){
                check_time=corrts;
                //! fill neutron here
                BELENHit* neuhitclone=new BELENHit;
                neuhit->Copy(*neuhitclone);
                if ((corrts-(Long64_t)ts)>=fNeuBetaoffset) flocalbetaS->AddNeutronForward(neuhitclone);
                else flocalbetaS->AddNeutronBackward(neuhitclone);
                ndelayedneutron++;
            }
            fhe3Map_it++;
        }
        if (ndelayedneutron>0) ncorrwneutron++;

        //!******************************ANCINARY detector****************

        /*
        //! Correlate imp with f11 right side only
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
        */

        //! Correlate imp with f11 left and right
        ts1 = (Long64_t)ts - (Long64_t)fBetaAncTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaAncTWup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fF11MapLR_it = fF11LRMap.lower_bound(ts1);
        while(fF11MapLR_it!=fF11LRMap.end()&&fF11MapLR_it->first<ts2){
            corrts = (Long64_t) fF11MapLR_it->first;
            correntry = fF11MapLR_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                ftrAnc->GetEvent(correntry);
                BELENHit* hit=new BELENHit;
                fanc->Copy(*hit);
                flocalbetaS->AddAnc(hit);
                ncorr++;
                break;
            }
            fF11MapLR_it++;
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

        //! Correlate with gamma
        ts1 = (Long64_t)ts - (Long64_t)fBetaGammaTWlow;
        ts2 = (Long64_t)ts + (Long64_t)fBetaGammaTWup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fcloverMap_it = fcloverMap.lower_bound(ts1);
        while(fcloverMap_it!=fcloverMap.end()&&fcloverMap_it->first<ts2){
            corrts = (Long64_t) fcloverMap_it->first;
            correntry = fcloverMap_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                ftrGamma->GetEvent(correntry);
                CloverHit* hit=new CloverHit;
                fclover->Copy(*hit);
                flocalbetaS->AddClover(hit);
                ncorr++;
                //break;
            }
            fcloverMap_it++;
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
        faidaImplantMapFull_it = faidaImplantMapFull.lower_bound(ts1);

        while(faidaImplantMapFull_it!=faidaImplantMapFull.end()&&faidaImplantMapFull_it->first<ts2){
            corrts = (Long64_t) faidaImplantMapFull_it->first;
            IonBetaMult* imp = faidaImplantMapFull_it->second;
            double impx= imp->GetHitPositionX();
            double impy= imp->GetHitPositionY();
            short impz= imp->GetHitPositionZ()+imp->GetDZ(); // corrected dZ

            //! stuffs for rejecting noise events associated with implant
            if (corrts!=check_time){// avoid multiple filling corrts!=check_time
                if (((Long64_t)ts-corrts)<fimpnoisefilter_dt&&((Long64_t)ts-corrts)>0){
                    if (((betax-impx>=-fimpnoisefilter_dxy)&&(betax-impx<=fimpnoisefilter_dxy)&&(betay-impy>=-fimpnoisefilter_dxy)&&(betay-impy<=fimpnoisefilter_dxy))&&betaz<=impz+fimpnoisefilter_dz){
                        flocalbetaS->SetDtIonAll((Double_t)((Long64_t)ts-corrts)/1e3);//unit us
                    }
                }
            }

            if (corrts!=check_time&&betaz==impz){// avoid multiple filling corrts!=check_time
                if (!((betax-impx>-fmaxnpixels)&&(betax-impx<fmaxnpixels)&&(betay-impy>-fmaxnpixels)&&(betay-impy<fmaxnpixels))){
                    faidaImplantMapFull_it++;
                    continue;
                }
                check_time=corrts;

                IonBetaMult* imparr=(IonBetaMult*)flocalimparray->ConstructedAt(nionbetacorr);

                imp->CopyWithBigRIPSOnly(*imparr);

                nionbetacorr++;
            }
            faidaImplantMapFull_it++;
        }


        //! fill beta data here
        flocalbetaS->CopyFromAIDA(beta);
        ftree->Fill();

        if (nionbetacorr>0) ncorrwbeta++;
        k++;
    }



    TSpectrum *s = new TSpectrum();
    ftotaltimepulser=(Double_t)(ftsendpulser-ftsbeginpulser)/1e9;

    //Double_t totaltime=(Double_t)(ftsendveto-ftsbeginveto)/1e9;
    Double_t totaltime=ftotaltimepulser;
    Double_t expectedcounts=totaltime*10;
    Double_t ncountsDtPuser=fhpulserall[140]->GetEntries();
    cout<<expectedcounts<<"/-/"<<ncountsDtPuser<<"/-/"<<fh1dtpulser->GetEntries()<<endl;

    //! new calculation of dead time
    for (Int_t i=0;i<141;i++){
        s->Search(fhpulser[i]);
        Double_t *xpeaks = s->GetPositionX();
        Double_t xp=0;
        if (s->GetNPeaks()>0) {
            xp = xpeaks[0];
            fhpulser[i]->Fit("gaus","RQ","",xp-100,xp+100);
            Double_t sigma=fhpulser[i]->GetFunction("gaus")->GetParameter(2);
            Double_t totalcounts=fhpulser[i]->Integral(fhpulser[i]->GetXaxis()->FindBin(xp-sigma*10),fhpulser[i]->GetXaxis()->FindBin(xp+sigma*10));
            Double_t totalcountsall=fhpulserall[i]->Integral(fhpulserall[i]->GetXaxis()->FindBin(xp-sigma*10),fhpulserall[i]->GetXaxis()->FindBin(xp+sigma*10));
            //cout<<i<<"-"<<totalcounts<<"-"<<totalcountsall<<endl;
            fh2deadtime->Fill(i,100-totalcounts/expectedcounts*100);
            fh1deadtime->Fill(100-totalcounts/expectedcounts*100);

            fh2deadtime2->Fill(i,100-totalcountsall/expectedcounts*100);
            fh1deadtime2->Fill(100-totalcountsall/expectedcounts*100);

            fh2deadtime3->Fill(i,100-totalcounts/ncountsDtPuser*100);
            fh1deadtime3->Fill(100-totalcounts/ncountsDtPuser*100);

            fh2deadtime4->Fill(i,100-totalcountsall/ncountsDtPuser*100);
            fh1deadtime4->Fill(100-totalcountsall/ncountsDtPuser*100);
        }
    }
    delete s;

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

void Merger::DoMergeClosePixel()
{

    for (Long64_t jentry=0;jentry<fnentriesMerged;jentry++){
        ftrMerged->GetEntry(jentry);
        //! ****** BETA VETO ***************
        //! downstream veto cut
        Int_t ndownstreamveto=0;
        for (Int_t j=0;j<flocalbetaS->GetAncMultipliticy();j++){
            if (flocalbetaS->GetAncHit(j)->GetMyPrecious()==4) ndownstreamveto++;
        }
        //! f11 cut
        Int_t nf11rveto=0;
        for (Int_t j=0;j<flocalbetaS->GetAncMultipliticy();j++){
            if (flocalbetaS->GetAncHit(j)->GetMyPrecious()==1||flocalbetaS->GetAncHit(j)->GetID()==1) nf11rveto++;
        }

        //! beta cut
        //if (flocalbetaS->GetDtIon()<=dtioncut||flocalbetaS->GetSumEXYRank()>sumexyrankcut||ndownstreamveto>0) continue;
        if (flocalbetaS->GetDtIonAll()>=0||flocalbetaS->GetSumEXYRank()>sumexyrankcut||ndownstreamveto>0) continue;

        if (!abs((int)((long long)flocalbetaS->GetXTimestamp()-(long long)flocalbetaS->GetYTimestamp())<xytdiffcut)) continue;// time cut
        //! select on multiplicity 1 events
        //if (!(flocalbetaS->GetXClusterMultiplicity()==1&&flocalbetaS->GetYClusterMultiplicity()==1)) continue;

        AIDASimpleStruct* hit=new AIDASimpleStruct;
        hit->SetEventNumber(flocalbetaS->GetEventNumber());
        hit->SetHitPosition(flocalbetaS->GetHitPositionX(),flocalbetaS->GetHitPositionY(),flocalbetaS->GetHitPositionZ());
        fcorrbetaMap.insert(make_pair(flocalbetaS->GetTimestamp(),make_pair(jentry,hit)));
    }
    cout<<"start correlation"<<endl;

    for (Long64_t jentry=0;jentry<fnentriesImp;jentry++){
        ftreeImplant->GetEvent(jentry);
        unsigned long long ts=flocalimp->GetTimestamp();
        double impx= flocalimp->GetHitPositionX();
        double impy= flocalimp->GetHitPositionY();
        short impz= flocalimp->GetHitPositionZ()+flocalimp->GetDZ();//! corrected dZ

        //! with beta
        Long64_t ts1 = ts - fIonBetaTWup;
        Long64_t ts2 = ts + fIonBetaTWlow;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fcorrbetaMap_it = fcorrbetaMap.lower_bound(ts1);

        unsigned int ievt = 0;
        unsigned int ievt_prev = 0;
        unsigned int betaclosepixelentry=0;
        unsigned long long betaclosepixelts=0;
        Double_t betaclosepixeldistance = fmaxnpixels*2;
        Double_t distance = 0;

        while(fcorrbetaMap_it!=fcorrbetaMap.end()&&fcorrbetaMap_it->first<ts2){
            corrts = fcorrbetaMap_it->first;
            correntry = fcorrbetaMap_it->second.first;
            AIDASimpleStruct* hit = fcorrbetaMap_it->second.second;
            double betax=hit->GetHitPositionX();
            double betay=hit->GetHitPositionY();
            short betaz=hit->GetHitPositionZ();

            if (corrts!=check_time&&betaz==impz){
                if (!((betax-impx>-deltaxcut)&&(betax-impx<deltaxcut)&&(betay-impy>-deltaycut)&&(betay-impy<deltaycut))){
                    fcorrbetaMap_it++;
                    continue;
                }
                ievt = hit->GetEventNumber();
                distance = sqrt((betax-impx)*(betax-impx)+(betay-impy)*(betay-impy));
                if (ievt==ievt_prev&&ievt_prev!=0){
                    if (distance<betaclosepixeldistance){
                        betaclosepixeldistance = distance;
                        betaclosepixelentry = correntry;
                        betaclosepixelts = corrts;
                    }
                }else{//! next evt
                    if (ievt!=ievt_prev){
                        ftrMerged->GetEntry(betaclosepixelentry);
                        //! fill events
                        ResetSimpleData();
                        decay.evt=flocalbetaS->GetEventNumber();
                        decay.ts=flocalbetaS->GetTimestamp();
                        decay.x=flocalbetaS->GetHitPositionX();
                        decay.y=flocalbetaS->GetHitPositionY();
                        decay.z=flocalbetaS->GetHitPositionZ();
                        decay.ex=flocalbetaS->GetXEnergy();
                        decay.ey=flocalbetaS->GetYEnergy();
                        decay.multx=flocalbetaS->GetXMultiplicity();
                        decay.multy=flocalbetaS->GetYMultiplicity();
                        decay.multz=flocalbetaS->GetZMultiplicity();

                        //! ****** NEUTRON ***************
                        //! neutron multipliciy counters and cut
                        Int_t nneufreal=0;
                        Int_t nneubreal=0;
                        for (Int_t j=0;j<flocalbetaS->GetNeutronForwardMultipliticy();j++){
                            if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetFinalVetoTime()<0){
                            //if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetF11Time()<0){
                                decay.neu_ch[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetID();
                                decay.neu_E[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetEnergy();
                                decay.neu_T[nneufreal]=((Double_t)((Long64_t)flocalbetaS->GetNeutronForwardHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;
                                decay.neu_x[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().X();
                                decay.neu_y[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().Y();
                                nneufreal++;
                            }
                        }

                        for (Int_t j=0;j<flocalbetaS->GetNeutronBackwardMultipliticy();j++){
                            if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetFinalVetoTime()<0){
                            //if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetF11Time()<0){
                                nneubreal++;
                            }
                        }
                        decay.neu_hit=nneufreal;
                        decay.neub_hit=nneubreal;

                        //!************ GAMMA *************
                        for (Int_t j=0;j<flocalbetaS->GetCloverMultipliticy();j++){
                            decay.gc_ch[j]=flocalbetaS->GetCloverHit(j)->GetCloverLeaf()+(flocalbetaS->GetCloverHit(j)->GetClover()-1)*4;
                            decay.gc_E[j]=flocalbetaS->GetCloverHit(j)->GetEnergy();
                            decay.gc_T[j]=((Double_t)((Long64_t)flocalbetaS->GetCloverHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;

                            //! gamma calibration goes here
                            //convert back to adc
                            decay.gc_E[j]=(decay.gc_E[j]-fcoffsetold[decay.gc_ch[j]-1])/fcgainold[decay.gc_ch[j]-1];
                            //apply new calibration
                            if (decay.gc_E[j]<fsep[decay.gc_ch[j]-1]){//low energy calibration
                                decay.gc_E[j]=flow_offset[decay.gc_ch[j]-1]+flow_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+flow_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
                            }else{//high energy calibration
                                decay.gc_E[j]=fhigh_offset[decay.gc_ch[j]-1]+fhigh_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+fhigh_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
                            }
                        }
                        decay.gc_hit=flocalbetaS->GetCloverMultipliticy();

                        //! implanation hit
                        decay.t=((Double_t)((Long64_t)flocalbetaS->GetTimestamp()-(Long64_t)flocalimp->GetTimestamp()))/1e9;
                        decay.ion_x=flocalimp->GetHitPositionX();
                        decay.ion_y=flocalimp->GetHitPositionY();
                        decay.ion_z=flocalimp->GetHitPositionZ();
                        decay.zet=flocalimp->GetBeamHit()->zet;

                        decay.aoq=flocalimp->GetBeamHit()->aoq;
                        decay.deltaxy=sqrt(decay.x*decay.ion_x+decay.y*decay.ion_y);
                        decay.ndecay=nionbetacorr;
                        decay.nbeta=nionbetacorr;
                        double zet=flocalimp->GetBeamHit()->zet;
                        double aoq=flocalimp->GetBeamHit()->aoq;
                        for (Int_t j=0;j<nri;j++){
                            if (!enablepid2[j]) continue;
                            if (cutg[j]->IsInside(aoq,zet)){
                                ftreeRI[j]->Fill();
                            }
                        }
                        ftree->Fill();
                        fh2->Fill((Long64_t)betaclosepixelts-(Long64_t)ts);
                    }
                    betaclosepixeldistance = distance;
                    betaclosepixelentry = correntry;
                    betaclosepixelts = corrts;
                }
                ievt_prev = ievt;
                fh1->Fill((Long64_t)corrts-(Long64_t)ts);
                //ftrMerged->GetEntry(correntry);
                ncorr++;
            }
            fcorrbetaMap_it++;
        }
        if (ncorr>0) {
            ftrMerged->GetEntry(betaclosepixelentry);
            //! fill events
            ResetSimpleData();
            decay.evt=flocalbetaS->GetEventNumber();
            decay.ts=flocalbetaS->GetTimestamp();
            decay.x=flocalbetaS->GetHitPositionX();
            decay.y=flocalbetaS->GetHitPositionY();
            decay.z=flocalbetaS->GetHitPositionZ();
            decay.ex=flocalbetaS->GetXEnergy();
            decay.ey=flocalbetaS->GetYEnergy();
            decay.multx=flocalbetaS->GetXMultiplicity();
            decay.multy=flocalbetaS->GetYMultiplicity();
            decay.multz=flocalbetaS->GetZMultiplicity();

            //! ****** NEUTRON ***************
            //! neutron multipliciy counters and cut
            Int_t nneufreal=0;
            Int_t nneubreal=0;

            for (Int_t j=0;j<flocalbetaS->GetNeutronForwardMultipliticy();j++){
                if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetFinalVetoTime()<0){
                //if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetF11Time()<0){
                    decay.neu_ch[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetID();
                    decay.neu_E[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetEnergy();
                    decay.neu_T[nneufreal]=((Double_t)((Long64_t)flocalbetaS->GetNeutronForwardHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;
                    decay.neu_x[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().X();
                    decay.neu_y[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().Y();
                    nneufreal++;
                }
            }

            for (Int_t j=0;j<flocalbetaS->GetNeutronBackwardMultipliticy();j++){
                if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetFinalVetoTime()<0){
                //if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetF11Time()<0){
                    nneubreal++;
                }
            }
            decay.neu_hit=nneufreal;
            decay.neub_hit=nneubreal;

            //!************ GAMMA *************
            for (Int_t j=0;j<flocalbetaS->GetCloverMultipliticy();j++){
                decay.gc_ch[j]=flocalbetaS->GetCloverHit(j)->GetCloverLeaf()+(flocalbetaS->GetCloverHit(j)->GetClover()-1)*4;
                decay.gc_E[j]=flocalbetaS->GetCloverHit(j)->GetEnergy();
                decay.gc_T[j]=((Double_t)((Long64_t)flocalbetaS->GetCloverHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;

                //! gamma calibration goes here
                //convert back to adc
                decay.gc_E[j]=(decay.gc_E[j]-fcoffsetold[decay.gc_ch[j]-1])/fcgainold[decay.gc_ch[j]-1];
                //apply new calibration
                if (decay.gc_E[j]<fsep[decay.gc_ch[j]-1]){//low energy calibration
                    decay.gc_E[j]=flow_offset[decay.gc_ch[j]-1]+flow_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+flow_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
                }else{//high energy calibration
                    decay.gc_E[j]=fhigh_offset[decay.gc_ch[j]-1]+fhigh_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+fhigh_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
                }
            }
            decay.gc_hit=flocalbetaS->GetCloverMultipliticy();

            //! implanation hit
            decay.t=((Double_t)((Long64_t)flocalbetaS->GetTimestamp()-(Long64_t)flocalimp->GetTimestamp()))/1e9;
            decay.ion_x=flocalimp->GetHitPositionX();
            decay.ion_y=flocalimp->GetHitPositionY();
            decay.ion_z=flocalimp->GetHitPositionZ();
            decay.zet=flocalimp->GetBeamHit()->zet;


            decay.aoq=flocalimp->GetBeamHit()->aoq;
            decay.deltaxy=sqrt(decay.x*decay.ion_x+decay.y*decay.ion_y);
            decay.ndecay=nionbetacorr;
            decay.nbeta=nionbetacorr;
            double zet=flocalimp->GetBeamHit()->zet;
            double aoq=flocalimp->GetBeamHit()->aoq;
            for (Int_t j=0;j<nri;j++){
                if (!enablepid2[j]) continue;
                if (cutg[j]->IsInside(aoq,zet)){
                    ftreeRI[j]->Fill();
                }
            }
            ftree->Fill();
            fh2->Fill((Long64_t)betaclosepixelts-(Long64_t)ts);
        }

    }

}

void Merger::DoSeparatePIDFinalTree()
{

    for (Long64_t jentry=0;jentry<fnentriesMerged;jentry++){
        ftrMerged->GetEntry(jentry);
        //! ****** BETA VETO ***************
        //! downstream veto cut
        Int_t ndownstreamveto=0;
        for (Int_t j=0;j<flocalbetaS->GetAncMultipliticy();j++){
            if (flocalbetaS->GetAncHit(j)->GetMyPrecious()==4) ndownstreamveto++;
        }
        //! f11 cut
        Int_t nf11rveto=0;
        for (Int_t j=0;j<flocalbetaS->GetAncMultipliticy();j++){
            if (flocalbetaS->GetAncHit(j)->GetMyPrecious()==1||flocalbetaS->GetAncHit(j)->GetID()==1) nf11rveto++;
        }

        //! beta cut        
        //if (flocalbetaS->GetDtIon()<=dtioncut||flocalbetaS->GetSumEXYRank()>sumexyrankcut||ndownstreamveto>0) continue;
        if (flocalbetaS->GetDtIonAll()>=0||flocalbetaS->GetSumEXYRank()>sumexyrankcut||ndownstreamveto>0) continue;
        if (!abs((int)((long long)flocalbetaS->GetXTimestamp()-(long long)flocalbetaS->GetYTimestamp())<xytdiffcut)) continue;// time cut
        //! select on multiplicity 1 events
        //if (!(flocalbetaS->GetXClusterMultiplicity()==1&&flocalbetaS->GetYClusterMultiplicity()==1)) continue;

        ResetSimpleData();
        decay.evt=flocalbetaS->GetEventNumber();
        decay.ts=flocalbetaS->GetTimestamp();
        decay.x=flocalbetaS->GetHitPositionX();
        decay.y=flocalbetaS->GetHitPositionY();
        decay.z=flocalbetaS->GetHitPositionZ();
        decay.ex=flocalbetaS->GetXEnergy();
        decay.ey=flocalbetaS->GetYEnergy();
        decay.multx=flocalbetaS->GetXMultiplicity();
        decay.multy=flocalbetaS->GetYMultiplicity();
        decay.multz=flocalbetaS->GetZMultiplicity();

        //! ****** NEUTRON ***************
        //! neutron multipliciy counters and cut
        Int_t nneufreal=0;
        Int_t nneubreal=0;


        for (Int_t j=0;j<flocalbetaS->GetNeutronForwardMultipliticy();j++){
            if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetFinalVetoTime()<0){
            //if (flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronForwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronForwardHit(j)->GetF11Time()<0){
                decay.neu_ch[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetID();
                decay.neu_E[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetEnergy();
                decay.neu_T[nneufreal]=((Double_t)((Long64_t)flocalbetaS->GetNeutronForwardHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;
                decay.neu_x[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().X();
                decay.neu_y[nneufreal]=flocalbetaS->GetNeutronForwardHit(j)->GetRndPosition().Y();
                nneufreal++;
            }
        }

        for (Int_t j=0;j<flocalbetaS->GetNeutronBackwardMultipliticy();j++){
            if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetFinalVetoTime()<0){
            //if (flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()>neutronecut[0]&&flocalbetaS->GetNeutronBackwardHit(j)->GetEnergy()<neutronecut[1]&&flocalbetaS->GetNeutronBackwardHit(j)->GetF11Time()<0){
                nneubreal++;
            }
        }
        decay.neu_hit=nneufreal;
        decay.neub_hit=nneubreal;

        //!************ GAMMA *************
        for (Int_t j=0;j<flocalbetaS->GetCloverMultipliticy();j++){
            decay.gc_ch[j]=flocalbetaS->GetCloverHit(j)->GetCloverLeaf()+(flocalbetaS->GetCloverHit(j)->GetClover()-1)*4;
            decay.gc_E[j]=flocalbetaS->GetCloverHit(j)->GetEnergy();
            decay.gc_T[j]=((Double_t)((Long64_t)flocalbetaS->GetCloverHit(j)->GetTimestamp()-(Long64_t)flocalbetaS->GetTimestamp()))/1e3;

            //! gamma calibration goes here
            //convert back to adc
            decay.gc_E[j]=(decay.gc_E[j]-fcoffsetold[decay.gc_ch[j]-1])/fcgainold[decay.gc_ch[j]-1];
            //apply new calibration
            if (decay.gc_E[j]<fsep[decay.gc_ch[j]-1]){//low energy calibration
                decay.gc_E[j]=flow_offset[decay.gc_ch[j]-1]+flow_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+flow_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
            }else{//high energy calibration
                decay.gc_E[j]=fhigh_offset[decay.gc_ch[j]-1]+fhigh_gain[decay.gc_ch[j]-1]*decay.gc_E[j]+fhigh_se[decay.gc_ch[j]-1]*decay.gc_E[j]*decay.gc_E[j];
            }
        }
        decay.gc_hit=flocalbetaS->GetCloverMultipliticy();

        //!************ Implantation vector *************
        Int_t nimpevt=flocalimparray->GetEntriesFast();
        for (Int_t j=0;j<nimpevt;j++){
            IonBetaMult* imp=(IonBetaMult*)flocalimparray->At(j);

            Double_t deltax=flocalbetaS->GetHitPositionX()-imp->GetHitPositionX();
            Double_t deltay=flocalbetaS->GetHitPositionX()-imp->GetHitPositionX();
            if(deltax*deltax>deltaxcut*deltaxcut||deltay*deltay>deltaycut*deltaycut) continue;

            decay.t=((Double_t)((Long64_t)flocalbetaS->GetTimestamp()-(Long64_t)imp->GetTimestamp()))/1e9;
            decay.ion_x=imp->GetHitPositionX();
            decay.ion_y=imp->GetHitPositionY();
            decay.ion_z=imp->GetHitPositionZ();
            decay.zet=imp->GetBeamHit()->zet;

            //! only for ahn experiment
            /*
            decay.zet=-9999;
            if (imp->GetBeamHit()->f11x>0) decay.zet=imp->GetBeamHit()->f11x;
            if (imp->GetBeamHit()->f11y>0) decay.zet=imp->GetBeamHit()->f11y;
            */

            decay.aoq=imp->GetBeamHit()->aoq;
            decay.deltaxy=sqrt(decay.x*decay.ion_x+decay.y*decay.ion_y);
            decay.ndecay=nionbetacorr;
            decay.nbeta=nionbetacorr;
            double zet=imp->GetBeamHit()->zet;
            double aoq=imp->GetBeamHit()->aoq;
            for (Int_t j=0;j<nri;j++){
                if (!enablepid2[j]) continue;
                if (cutg[j]->IsInside(aoq,zet)){
                    ftreeRI[j]->Fill();
                }
            }
            ftree->Fill();
        }
    }
    cout<<"Making implant tree ..."<<endl;

    Int_t nimplant[nri];
    for (Int_t i=0;i<nri;i++) nimplant[i]=0;

    Long64_t tsstart=0;
    Long64_t tsstop=0;

    //! making separated implant tree
    for (Long64_t jentry=0;jentry<fnentriesImp;jentry++){
        ftreeImplant->GetEvent(jentry);
        double zet=flocalimp->GetBeamHit()->zet;
        double aoq=flocalimp->GetBeamHit()->aoq;
        if (zet>0&&aoq>0){
            if (tsstart==0) tsstart=flocalimp->GetTimestamp();
            tsstop=flocalimp->GetTimestamp();
        }
        for (Int_t j=0;j<nri;j++){
            if (!enablepid2[j]) continue;
            if (cutg[j]->IsInside(aoq,zet)){                
                ftreeimplantRI[j]->Fill();
                if (flocalimp->GetHitPositionZ()<5) nimplant[j]++;
            }
        }
        ftreeimplantAll->Fill();
    }

    std::ofstream ofs("out.txt",ios::app);
    ofs<<finputMerged<<"\t";
    for (Int_t i=0;i<nri;i++){
        if (!enablepid2[i]) continue;
        ofs<<nimplant[i]<<"\t";
    }
    ofs<<tsstop-tsstart<<endl;

}
