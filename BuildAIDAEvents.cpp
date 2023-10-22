#include "BuildAIDAEvents.h"
BuildAIDAEvents::BuildAIDAEvents()
{
    faidatranst = 20000;
    fwindowBeta = 2500;
    fwindowIon = 2500;
    fwindowDisc = 5000;
    fflag_filldata = false;
    ftemp = 0;
    fflag_pulser_in_stream = false;
    fflag_corrscaler_in_stream = false;

    for (Int_t i=0;i<NumDSSD;i++){
        fsumexcut[i] = 0;
        fsumeycut[i] = 0;
    }
    fcorrcut = -1;
    fisgzstream = false;
    aidaunpkg = new AIDAUnpacker;


    flocalaidaBETA = new AIDA;
    flocalaidaION = new AIDA;
    flocalaidaCORR = new AIDA;
}

BuildAIDAEvents::~BuildAIDAEvents()
{
    //! error on this
    //delete aidaunpkg;
    delete flocalaidaBETA;
    delete flocalaidaION;
    delete flocalaidaCORR;

}

void BuildAIDAEvents::Init(char* aidafile)
{
    flastevent = -1;

    fflag_addFirstIonHit=false;
    fflag_addFirstBetaHit=false;
    fisbeta = true;

    fADIonEntry = 0;
    fADBetaEntry = 0;
    fADPulserEntry = 0;
    fADts = 0;

    //! for writing
    faida = new AIDA;

    //! initilize input aida
    if (fisgzstream) aidaunpkg->SetGzStream();
    aidaunpkg->Init(aidafile);
    aidaunpkg->read_mapping(fmappingfile);
    if (fthresholdfile!=NULL) aidaunpkg->read_threshold_table(fthresholdfile);
    aidaunpkg->SetVerbose(fverbose);
    cout<<"Trying to get first Sync...."<<endl;
    aidaunpkg->GetFirstSync();
    cout<<"\n\n****** Start reading AIDA ******\n\n"<<endl;

    //! Get number of entry
    //!TEMPORARY
    fADNblock = aidaunpkg->GetNBlock();
    cout << fADNblock << " blocks in AIDA file" << endl;

    //! reset last time stamp
    flastADIonts = 0;
    fcurrentionts = 0;
    flastADBetats = 0;
    fcurrentbetats = 0;
    for (Int_t i=0;i<NumFee;i++) {
        for (Int_t j=0;j<NumChFee;j++){
            flastfastts[i][j] = 0;
            flastfasttsEXT[i][j] = 0;
        }
    }
    fflag_trans = false;
}

void BuildAIDAEvents::BookTree(TTree *treeIon,TTree *treeBeta,TTree *treePulser,Int_t bufsize)
{
    //! initilize output
    fmtrION = treeIon;
    fmtrION->Branch("adentry",&fADIonEntry,bufsize); //320000
    fmtrION->Branch("adTS",&fADtsION,bufsize);
    fmtrION->Branch("aida",&flocalaidaION,bufsize);
    fmtrION->BranchRef();


    fmtrBETA = treeBeta;
    fmtrBETA->Branch("adentry",&fADBetaEntry,bufsize);
    fmtrBETA->Branch("adTS",&fADtsBETA,bufsize);
    fmtrBETA->Branch("aida",&flocalaidaBETA,bufsize);
    fmtrBETA->BranchRef();

    fmtrPULSER = treePulser;
    fmtrPULSER->Branch("adentry",&fADPulserEntry,bufsize);
    fmtrPULSER->Branch("adTS",&fADtsBETA,bufsize);
    fmtrPULSER->Branch("aida",&flocalaidaBETA,bufsize);
    fmtrPULSER->BranchRef();

    fflag_filldata=true;
}

void BuildAIDAEvents::ReadCalibTable()
{

    //clean up
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_cal[i][j][0]=0.;
            dssd_cal[i][j][1]=1.;
        }
    }

    ifstream inpf(fcalibfile);
    if (inpf.fail()){
        cout<<"No Calibration table is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table: "<<fcalibfile<<endl;
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

//! Add AIDA ION hits
void BuildAIDAEvents::AddAIDAIonHits(rawaida_info aidaraw){
    AIDAHit* hit = new AIDAHit;
    hit->SetADC(aidaraw.adcData);
    hit->SetEnergy((double)aidaraw.adcData);
    hit->SetTimestamp(aidaraw.extTimestamp);
    //hit->SetTimestamp(aidaraw.timestamp);
    hit->SetID(aidaraw.stripNo + aidaraw.dssdNo*NumStrXY);
    hit->SetXY(aidaraw.stripNo);
    hit->SetZ(aidaraw.dssdNo);
    hit->SetFEE(aidaraw.feeNo);
    hit->SetFEEChannel(aidaraw.chNo);
    flocalaidaION->AddHit(hit);
}

//! Add AIDA Beta hitsflocalbeta
void BuildAIDAEvents::AddAIDABetaHits(rawaida_info aidaraw){
    AIDAHit* hit = new AIDAHit;
    hit->SetADC(aidaraw.adcData);
    hit->SetEnergy((double)aidaraw.adcData*dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][0]);
    hit->SetTimestamp(aidaraw.extTimestamp);
    //hit->SetTimestamp(aidaraw.timestamp);
    if ((flastfastts[aidaraw.feeNo][aidaraw.chNo] + fwindowDisc) > aidaraw.timestamp)
        hit->SetFastTimestamp(flastfasttsEXT[aidaraw.feeNo][aidaraw.chNo]);
    hit->SetID(aidaraw.stripNo + aidaraw.dssdNo*NumStrXY);
    hit->SetXY(aidaraw.stripNo);
    hit->SetZ(aidaraw.dssdNo);
    hit->SetFEE(aidaraw.feeNo);
    hit->SetFEEChannel(aidaraw.chNo);
    flocalaidaBETA->AddHit(hit);
}

//!close event ION
bool BuildAIDAEvents::CloseIonEvent()
{
    if(fverbose>0)
      cout << __PRETTY_FUNCTION__ << endl;
    //cout<<fADIonEntry<<"-"<<flocalaidaION->GetMult()<<endl;
    //fADIonEntry++;
    if (flocalaidaION->IonGetPos()) {
        flocalaidaION->SetTimestamp(flocalaidaION->GetHit(0)->GetTimestamp());
        if (fflag_filldata) fmtrION->Fill();
        fADIonEntry++;
        return true;
    }
    return false;
}

//!close event BETA
bool BuildAIDAEvents::CloseBetaEvent()
{
    if(fverbose>0)
      cout << __PRETTY_FUNCTION__ << endl;

    unsigned short* hitx= flocalaidaBETA->GetMultXs();
    unsigned short* hity= flocalaidaBETA->GetMultYs();
    unsigned short mult=flocalaidaBETA->GetMult();

    if (mult<64&&!fflag_trans&&(hitx[0]+hity[0]<8)&&(hitx[1]+hity[1]<8)&&(hitx[2]+hity[2]<8)&&(hitx[3]+hity[3]<8)&&(hitx[4]+hity[4]<8)&&(hitx[5]+hity[5]<8))
    {
        if (flocalaidaBETA->BetaGetPosNew(fcorrcut,fsumexcut,fsumeycut)){
            flocalaidaBETA->SetTimestamp(flocalaidaBETA->GetHit(0)->GetTimestamp());
            if (fflag_filldata) fmtrBETA->Fill();
            fADBetaEntry++;
            return true;
        }else {
            return false;
        }
    }
    if (mult>=64){
        flocalaidaBETA->SetTimestamp(fADtsBETA);
        if (fflag_filldata) fmtrPULSER->Fill();
        fADPulserEntry++;
        if (fflag_pulser_in_stream) return true;
        else return false;
        //!will not give to the "nextevent" scheme!
    }
    return false;
}

//!AIDA event builder also
bool BuildAIDAEvents::GetNextEvent(){
    //! read from aida
    bool flag_stop=false;

    if (fflag_addFirstIonHit) {
        flocalaidaION->Clear();
        AddAIDAIonHits(aidaraw);
    }

    if (fflag_addFirstBetaHit) {
        flocalaidaBETA->Clear();
        AddAIDABetaHits(aidaraw);
    }

    fflag_addFirstIonHit=false;
    fflag_addFirstBetaHit=false;
    while(!flag_stop){
        if (!aidaunpkg->GetNextHit()) return false;

        //! handle correlation scaler

        if (fflag_corrscaler_in_stream&&aidaraw.rangeType==-1){
            AIDAHit* hit = new AIDAHit;
            //if we dont set all things -> memory leak
            hit->SetADC(0);
            hit->SetEnergy(0.);
            hit->SetTimestamp(aidaraw.extTimestamp*tm_stp_scaler_ratio);
            hit->SetID(9999);
            hit->SetXY(0);
            hit->SetZ(0);
            hit->SetFEE(0);
            hit->SetFEEChannel(0);
            hit->SetRange(0);
            flocalaidaCORR->AddHit(hit);
        }


        aidaraw= aidaunpkg->GetAIDAraw();
        fADcurblock = aidaunpkg->GetCurrentBlock();
        //!ION events
        if (aidaraw.infoCode==0&&aidaraw.rangeType==1){
            flastADIonts = aidaraw.timestamp;
            if (flastADIonts <= fcurrentionts + fwindowIon){
                AddAIDAIonHits(aidaraw);
            }else{ //close event
                flag_stop=CloseIonEvent();
                //! Set eariest hit time!
                fcurrentionts = flastADIonts;
                fADtsION=aidaraw.extTimestamp;
                if (!flag_stop) {
                    flocalaidaION->Clear();
                    AddAIDAIonHits(aidaraw);
                }else{
                    fflag_addFirstIonHit=true;
                    fisbeta = false;
                }
            }
        }

        //! BETA events
        if (aidaraw.infoCode==0&&aidaraw.rangeType==0){
            flastADBetats = aidaraw.timestamp;
            if (flastADBetats <= fcurrentbetats + fwindowBeta){//within the window
                if (!fflag_trans) AddAIDABetaHits(aidaraw);
            }else{ //close event
                flag_stop=CloseBetaEvent();
                //!Check if previous ION entry is not within the window
                if ((flastADIonts+faidatranst)>flastADBetats) fflag_trans = true;
                else fflag_trans = false;
                fcurrentbetats = flastADBetats;
                //! Set eariest hit time!
                fADtsBETA=aidaraw.extTimestamp;
                if (!flag_stop) {
                    flocalaidaBETA->Clear();
                    AddAIDABetaHits(aidaraw);
                }else{
                    fflag_addFirstBetaHit=true;
                    fisbeta = true;
                }
            }
        }
        //! handle fast discriminator time stamp
        if (aidaraw.infoCode==6){
            int feen = aidaraw.feeNo;
            int pattern = 0xFFFF & aidaraw.adcData;
            int asicsn = ((aidaraw.adcData & 0xF0000) >> 16) -1; //(from 0 to 3)
            //cout<<"feeno"<<std::dec<<feen<<" adcdata 0x"<<std::hex<<(aidaraw.adcData&0xFFFFF)<<" asics no"<<std::dec<<asicsn<<std::hex<<" pattern 0x"<<pattern <<endl;
            for (int i=0;i<16;i++){
                if (((pattern >> i) &0x1) == 1) {
                    int chn = (asicsn*16)+i;
                    if (chn<NumChFee){
                        flastfastts[feen][chn] = aidaraw.timestamp;
                        flastfasttsEXT[feen][chn] = aidaraw.extTimestamp * tm_stp_scaler_ratio;
                        //cout<<std::dec<<feen<<"-"<<chn<<endl;
                    }else{
                        cout<<"Somethings wrong with the fast discriminator data!"<<endl;
                    }

                }
            }
        }

        //! handle ion event within beta window
        if (aidaraw.infoCode==0&&(aidaraw.rangeType==1||(aidaraw.rangeType==0&&aidaraw.adcData>25000))){
            fflag_trans=true;
        }
    }
    return true;
}

bool BuildAIDAEvents::check_channel_mask(Int_t dssdNo, Int_t stripNo, Int_t pattern,Int_t fee_num)
{
    Int_t check_feeNumber=aidaunpkg->DSSDtoFee[dssdNo][stripNo];
    Int_t check_chNumber=aidaunpkg->DSSDtoCh[dssdNo][stripNo];
    Int_t check_asics=check_chNumber/16;

    if (check_chNumber>64) {
        cout<<"something wrong in the mapping!"<<dssdNo<<"-"<<stripNo<<endl;
        return false;
    }
    check_chNumber=check_chNumber%16;
    Int_t asics_num=pattern&0xF0000>>16;
    if(check_feeNumber==fee_num&&asics_num==check_asics&&((pattern>>check_chNumber)&0x1==1))
        return true;
    else return false;
}




