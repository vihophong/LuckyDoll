#include "BelenReader.h"

void genRndCircle(Double_t &x,Double_t &y,Double_t a,Double_t b,Double_t xpos,Double_t ypos,Double_t R){
    if (b<a){
        Double_t temp=a;
        a=b;
        b=temp;
    }
    x=xpos+b*R*TMath::Cos(2*TMath::Pi()*a/b);
    y=ypos+b*R*TMath::Sin(2*TMath::Pi()*a/b);
}

BelenReader::BelenReader():rr()
{

    for (Int_t i=0;i<MaxID;i++){
        fHe3Id2posX[i]=0;
        fHe3Id2posY[i]=0;
        fHe3Id2posZ[i]=0;
        fHe3Id2diameter[i]=0;
    }

    for (Int_t i=0;i<MaxID;i++){
        fCrystalId2posX[i]=0;
        fCrystalId2posY[i]=0;
        fCrystalId2posZ[i]=0;
    }


    ftreedataNeuron = NULL;
    ftreedataGamma = NULL;
    ftreedataAnc = NULL;
    fflag_filldata = false;
}

BelenReader::~BelenReader()
{
    delete flocalNeutron;
    delete flocalGamma;
    delete flocalAnc;
}


void BelenReader::Init(char* belenfile){
    fBLAncEntry = 0;
    fBLGamEntry = 0;
    fBLAncEntry = 0;

    flocalNeutron = new BELENHit;
    flocalGamma = new CloverHit;
    flocalAnc = new BELENHit;

    finfile = new TFile(belenfile);
    ftree = (TTree*) finfile->Get("BelenTree");

    fnentries = ftree->GetEntries();
    cout<<"There are "<<fnentries<<" entries in Belen: "<< belenfile<<endl;

    //! brach tree
    ftree->SetBranchAddress("Neutrons",&ftreedataNeuron);
    ftree->SetBranchAddress("Gamma",&ftreedataGamma);
    ftree->SetBranchAddress("Ancillary",&ftreedataAnc);

    fcurentry = 0;
    fBLNeuEntry = 0;
    fBLGamEntry = 0;
    fBLAncEntry = 0;
    GetMapping();
    //if (!GetNextEvent()) exit(1);
}

void BelenReader::GetMapping(){
    std::ifstream inpf(fmappingfile);
    if (inpf.fail()){
        cout<<"No BELEN Mapping file is given"<<endl;
        return;
    }
    cout<<"Start reading BELEN Mapping file: "<<fmappingfile<<endl;

    Int_t He3id,daqId;
    Double_t x,y,z;
    Double_t d;
    Int_t mm=0;

    while (inpf.good()){
        inpf>>He3id>>daqId>>d>>x>>y>>z;
        fHe3Id2posX[daqId]=x;
        fHe3Id2posY[daqId]=y;
        fHe3Id2posZ[daqId]=z;
        fHe3Id2diameter[daqId]=d;
        //cout<<He3id<<"-"<<daqId<<"-"<<d<<"-"<<x<<"-"<<y<<"-"<<z<<endl;
        mm++;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    inpf.close();
}

void BelenReader::BookTree(TTree* treeNeutron, TTree *treeGamma, TTree *treeAnc, Int_t bufsize){
    //! initilize output
    fmtrNeutron = treeNeutron;
    fmtrNeutron->Branch("blentry",&fBLNeuEntry,bufsize); //320000
    fmtrNeutron->Branch("blTS",&fBLtsNeutron,bufsize);
    fmtrNeutron->Branch("belen",&flocalNeutron,bufsize);
    fmtrNeutron->BranchRef();


    fmtrGamma = treeGamma;
    fmtrGamma->Branch("blentry",&fBLGamEntry,bufsize);
    fmtrGamma->Branch("blTS",&fBLtsGamma,bufsize);
    fmtrGamma->Branch("belen",&flocalGamma,bufsize);
    fmtrGamma->BranchRef();

    fmtrAnc = treeAnc;
    fmtrAnc->Branch("blentry",&fBLAncEntry,bufsize);
    fmtrAnc->Branch("blTS",&fBLtsAnc,bufsize);
    fmtrAnc->Branch("belen",&flocalAnc,bufsize);
    fmtrAnc->BranchRef();

    fflag_filldata=true;
}


bool BelenReader::GetNextEvent(){
    ftree->GetEntry(fcurentry);

    fE = ftreedataNeuron->E + ftreedataGamma->E + ftreedataAnc->E;
    fT = ftreedataNeuron->T + ftreedataGamma->T + ftreedataAnc->T;
    fId = ftreedataNeuron->Id + ftreedataGamma->Id + ftreedataAnc->Id;
    ftype = ftreedataNeuron->type + ftreedataGamma->type + ftreedataAnc->type;
    fIndex1 = ftreedataNeuron->Index1 + ftreedataGamma->Index1 + ftreedataAnc->Index1;
    fIndex2 = ftreedataNeuron->Index2 + ftreedataGamma->Index2 + ftreedataAnc->Index2;
    fInfoFlag = ftreedataNeuron->InfoFlag + ftreedataGamma->InfoFlag + ftreedataAnc->InfoFlag;
    fName = ftreedataNeuron->Name + ftreedataGamma->Name + ftreedataAnc->Name;

    fcurentry++;

    //! fill data if nessacry

    if (ftype==1){
        flocalNeutron->SetEnergy(fE);
        flocalNeutron->SetTimestamp(fT);
        flocalNeutron->SetID(fId);
        PertubateHe3(fId);
        flocalNeutron->SetPos(fposX,fposY,fposZ);

        if (fflag_filldata) fmtrNeutron->Fill();
        fBLNeuEntry++;
    }else if (ftype==2){
        flocalGamma->SetEnergy(fE);
        flocalGamma->SetTimestamp(fT);
        flocalGamma->SetID(fId);
        PertubateClover(fId);
        flocalGamma->SetPos(fposX,fposY,fposZ);

        if (fflag_filldata) fmtrGamma->Fill();
        fBLGamEntry++;
    }
    else if (ftype==3){
        flocalAnc->SetEnergy(fE);
        flocalAnc->SetTimestamp(fT);
        flocalAnc->SetID(fId);
        if (fflag_filldata) fmtrAnc->Fill();
        fBLAncEntry++;
    }

    if (fcurentry>fnentries) return false;
    return true;
}

bool BelenReader::GetNextNeutronEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=1){
        if (!GetNextEvent()) return false;
    }
    fBLNeuEntry++;
    return true;
}

bool BelenReader::GetNextGammaEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=2){
        if (!GetNextEvent()) return false;
    }
    fBLGamEntry++;
    return true;
}

bool BelenReader::GetNextAncEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=3){
        if (!GetNextEvent()) return false;
    }
    fBLAncEntry++;
    return true;
}

void BelenReader::PertubateHe3(UShort_t He3Id){
    //!there is nothing here for the moment!
    fposX = fHe3Id2posX[He3Id];
    fposY = fHe3Id2posY[He3Id];
    fposZ = fHe3Id2posZ[He3Id];
    //! pertubating
    Double_t a,b,x,y,r;
    r=fHe3Id2diameter[He3Id]/2;
    a=rr.Rndm();
    b=rr.Rndm();
    genRndCircle(x,y,a,b,fposX,fposY,r);
    fposX = x;
    fposY = y;
}

void BelenReader::PertubateClover(UShort_t CrystalId){
    //!there is nothing here for the moment!
    fposX = fCrystalId2posX[CrystalId];
    fposY = fCrystalId2posZ[CrystalId];
    fposZ = fCrystalId2posZ[CrystalId];
    //! pertubating
}


