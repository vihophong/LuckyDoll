#include "BuildAIDAEvents.h"
BuildEvents::BuildEvents()
{
    faidatranst = 20000;
    fwindowBeta = 2500;
    fwindowIon = 2500;
}


BuildEvents::~BuildEvents()
{

}

void BuildEvents::Init(TTree* brtr, TTree* bltr, char* aidafile)
{

    flastevent = -1;

    fBRtr = brtr;
    fBLtr = bltr;

    fBRentry = 0;
    fBLentry = 0;
    fADIonEntry = 0;
    fADBetaEntry = 0;

    fBRts = 0;
    fADts = 0;
    fBLts = 0;

    //! for reading from tree and aida stream
    fBRbeam = NULL;
    fBReventinfo = NULL;
    fBLCLeventinfo = NULL;
    fBLCLhits = NULL;

    flocalaidaBETA = new AIDA;
    flocalaidaION = new AIDA;
    //! for writing
    faida = new AIDA;
    fbelen = new BELEN;
    fclover = new Clover;



    //! initilize input aida
    aidaunpkg = new AIDAUnpacker;
    aidaunpkg->Init(aidafile);
    aidaunpkg->read_mapping(fmappingfile);
    if (fthresholdfile!=NULL) aidaunpkg->read_threshold_table(fthresholdfile);
    aidaunpkg->SetVerbose(fverbose);
    cout<<"Trying to get first Sync...."<<endl;
    aidaunpkg->GetFirstSync();
    cout<<"\n\n****** Start reading AIDA ******\n\n"<<endl;

    //! initilize input bigrips
    fBRtr->SetBranchAddress("EventInfo",&fBReventinfo);
    fBRtr->SetBranchAddress("BigRIPSBeam",&fBRbeam);
    //! initilize input belen and clover
    //!
    //!

    //! Get number of entry
    //!TEMPORARY
    fADNblock = aidaunpkg->GetNBlock();
    fBRentries = fBRtr->GetEntries();

    cout << fBRentries << " entries in BigRIPS tree" << endl;
    cout << fADNblock << " blocks in AIDA file" << endl;
    //! Check first and last events for data integrity
    if(fverbose>0){
      Int_t status = fBRtr->GetEvent(0);
      if(status<0)
      cout << "first BigRIPS entry faulty!" << endl;
      //!newly added
      cout << "first BigRIPS timestamp: " << ((TArtEventInfo*) fBReventinfo->At(0))->GetTimeStamp() << endl;
      status = fBRtr->GetEvent(fBRentries-1);
      if(status<0)
        cout << "last BigRIPS entry faulty!" << endl;
      cout << "last BigRIPS timestamp: " << ((TArtEventInfo*) fBReventinfo->At(0))->GetTimeStamp() << endl;;
    }


    //! initilize output
    fmtrION = new TTree("ION","merged tree for ION");
    //fmtrION->Branch("brentry",&fBRentry,320000);
    //fmtrION->Branch("brTS",&fBRts,320000);
    //fmtrION->Branch("BRBeam",&fBRbeam,320000); //same as input
    fmtrION->Branch("adentry",&fADIonEntry,320000);
    fmtrION->Branch("adTS",&fADtsION,320000);
    fmtrION->Branch("aida",&flocalaidaION,320000);
    //fmtrION->Branch("blentry",&fBLentry,320000);
    //fmtrION->Branch("blTS",&fBLts,320000);
    //fmtrION->Branch("belen",&fbelen,320000);
    fmtrION->BranchRef();

    fmtrBETA = new TTree("BETA","merged tree for BETA");
    fmtrBETA->Branch("adentry",&fADBetaEntry,320000);
    fmtrBETA->Branch("adTS",&fADtsBETA,320000);
    fmtrBETA->Branch("aida",&flocalaidaBETA,320000);
    //fmtrBETA->Branch("blentry",&fBLentry,320000);
    //fmtrBETA->Branch("blTS",&fBLts,320000);
    //fmtrBETA->Branch("belen",&fbelen,320000);
    fmtrBETA->BranchRef();

    //! reset last time stamp
    flastBLts = 0;
    flastBRts = 0;
    flastADIonts = 0;
    fcurrentionts = 0;
    flastADBetats = 0;
    fcurrentbetats = 0;
    fflag_trans = false;
    //! clear detector
    fdetectors.clear();
}

void BuildEvents::ReadCalibTable()
{

    //clean up
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_cal[i][j][0]=0.;
            dssd_cal[i][j][1]=1.;
        }
    }

    cout<<"Reading calibration file: "<<fcalibfile<<endl;
    ifstream inpf(fcalibfile);
    if (inpf.fail()){
        cout<<"No Calibration table is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table"<<fcalibfile<<endl;
    Int_t dssd_index,strip_index;
    Double_t cal1,cal2;
    Int_t mm=0;

    while (inpf.good()){
    //for (Int_t i=0;i<100;i++){
        inpf>>dssd_index>>strip_index>>cal1>>cal2;
        dssd_cal[dssd_index][strip_index][0]=cal1;
        dssd_cal[dssd_index][strip_index][1]=cal2;
        //cout<<dssd_thr[dssd_index][strip_index]<<endl;
        mm++;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    inpf.close();


}

bool BuildEvents::ReadBigRIPS()
{
    if(fverbose>1)
      cout << __PRETTY_FUNCTION__ << endl;

    return true;
}

bool BuildEvents::ReadBELEN()
{
    if(fverbose>1)
      cout << __PRETTY_FUNCTION__ << endl;

    return true;
}

bool BuildEvents::ReadClover()
{
    if(fverbose>1)
      cout << __PRETTY_FUNCTION__ << endl;

    return true;
}

//! Add AIDA ION hits
void BuildEvents::AddAIDAIonHits(rawaida_info aidaraw){
    AIDAHit* hit = new AIDAHit;
    hit->SetADC(aidaraw.adcData);
    hit->SetEnergy(aidaraw.adcData * dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][0]);
    hit->SetTimestamp(aidaraw.extTimestamp);
    hit->SetID(aidaraw.dssdNo * aidaraw.stripNo);
    hit->SetXY(aidaraw.stripNo);
    hit->SetZ(aidaraw.dssdNo);
    flocalaidaION->AddHit(hit);
}

//! Add AIDA Beta hits
void BuildEvents::AddAIDABetaHits(rawaida_info aidaraw){
    AIDAHit* hit = new AIDAHit;
    hit->SetADC(aidaraw.adcData);
    hit->SetEnergy(aidaraw.adcData);
    hit->SetTimestamp(aidaraw.extTimestamp);
    hit->SetID(aidaraw.dssdNo * aidaraw.stripNo);
    hit->SetXY(aidaraw.stripNo);
    hit->SetZ(aidaraw.dssdNo);
    flocalaidaBETA->AddHit(hit);
}

//!close event ION
void BuildEvents::CloseIonEvent()
{
    if(fverbose>0)
      cout << __PRETTY_FUNCTION__ << endl;
    flocalaidaION->SetTimestamp(fADtsION);
    fmtrION->Fill();
    fADIonEntry++;
}

//!close event BETA
void BuildEvents::CloseBetaEvent()
{
    if(fverbose>0)
      cout << __PRETTY_FUNCTION__ << endl;

    unsigned short* hitx= flocalaidaBETA->GetMultX();
    unsigned short* hity= flocalaidaBETA->GetMultY();
    unsigned short mult=flocalaidaBETA->GetMult();

    if (mult<=64&&!fflag_trans&&(hitx[0]+hity[0]<8)&&(hitx[1]+hity[1]<8)&&(hitx[2]+hity[2]<8)&&(hitx[3]+hity[3]<8)&&(hitx[4]+hity[4]<8)&&(hitx[5]+hity[5]<8))
    {
        flocalaidaBETA->SetTimestamp(fADtsBETA);
        fmtrBETA->Fill();
        fADBetaEntry++;
    }
}

//!AIDA event builder also
bool BuildEvents::ReadAIDA()
{
    if(fverbose>1)
      cout << __PRETTY_FUNCTION__ << endl;

    //! read from aida
    if (!aidaunpkg->GetNextHit()) return false;
    fADcurblock = aidaunpkg->GetCurrentBlock();
    rawaida_info aidaraw= aidaunpkg->GetAIDAraw();

    //!ION events
    if (aidaraw.infoCode==0&&aidaraw.rangeType==1){
        flastADIonts = aidaraw.timestamp;
        if (flastADIonts <= fcurrentionts + fwindowIon){
            AddAIDAIonHits(aidaraw);
        }else{ //close event
            CloseIonEvent();
            flocalaidaION->Clear();
            fcurrentionts = flastADIonts;
            fADtsION=aidaraw.extTimestamp;
            //! Set eariest hit time!
            AddAIDAIonHits(aidaraw);
        }
    }


    //! BETA events
    if (aidaraw.infoCode==0&&aidaraw.rangeType==0){
        flastADBetats = aidaraw.timestamp;
        if (flastADBetats <= fcurrentbetats + fwindowBeta){
            if (!fflag_trans) AddAIDABetaHits(aidaraw);
        }else{ //close event
            CloseBetaEvent();
            //!Check if previous ION entry is not within the window

            if ((flastADIonts+faidatranst)>flastADBetats) fflag_trans = true;
            else fflag_trans = false;

            flocalaidaBETA->Clear();
            fcurrentbetats = flastADBetats;
            fADtsBETA=aidaraw.extTimestamp;
            //! Set eariest hit time!
            AddAIDABetaHits(aidaraw);
        }
    }

    if (aidaraw.rangeType==1||(aidaraw.rangeType==0&&aidaraw.adcData>25000)){
        fflag_trans=true;
    }

    return true;
}
bool BuildEvents::GetNextEvent(){
    //! read from aida
    bool flag_stop=false;
    while(!flag_stop){
        if (!aidaunpkg->GetNextHit()) return false;
        rawaida_info aidaraw= aidaunpkg->GetAIDAraw();
        fADcurblock = aidaunpkg->GetCurrentBlock();
        //!ION events
        if (aidaraw.infoCode==0&&aidaraw.rangeType==1){
            flastADIonts = aidaraw.timestamp;
            if (flastADIonts <= fcurrentionts + fwindowIon){
                AddAIDAIonHits(aidaraw);
            }else{ //close event
                CloseIonEvent();
                flocalaidaION->Clear();
                fcurrentionts = flastADIonts;
                fADtsION=aidaraw.extTimestamp;
                //! Set eariest hit time!
                AddAIDAIonHits(aidaraw);
                flag_stop=true;
            }
        }


        //! BETA events
        if (aidaraw.infoCode==0&&aidaraw.rangeType==0){
            flastADBetats = aidaraw.timestamp;
            if (flastADBetats <= fcurrentbetats + fwindowBeta){
                if (!fflag_trans) AddAIDABetaHits(aidaraw);
            }else{ //close event
                CloseBetaEvent();
                //!Check if previous ION entry is not within the window

                if ((flastADIonts+faidatranst)>flastADBetats) fflag_trans = true;
                else fflag_trans = false;

                flocalaidaBETA->Clear();
                fcurrentbetats = flastADBetats;
                fADtsBETA=aidaraw.extTimestamp;
                //! Set eariest hit time!
                AddAIDABetaHits(aidaraw);
                flag_stop=true;
            }
        }

        if (aidaraw.rangeType==1||(aidaraw.rangeType==0&&aidaraw.adcData>25000)){
            fflag_trans=true;
        }
    }
    return true;
}

bool BuildEvents::Merge()
{
    if (!ReadAIDA()) return false;
    //!BETA events
    return true;
}



