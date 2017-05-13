#include "BuildAIDAEventsNew.h"
BuildAIDAEvents::BuildAIDAEvents()
{
    fwindowHits = 200;
    flag_threhold = false;

    //faidatranst = 20000;
    //fwindowBeta = 2500;
    //fwindowIon = 2500;
    fwindowDisc = 5000;
    fflag_filldata = false;
    ftemp = 0;
    fflag_pulser_in_stream = false;
    fflag_corrscaler_in_stream = false;

    for (Int_t i=0;i<NumDSSD;i++){
        fsumexcut[i] = -1000.;
        fsumeycut[i] = -1000.;
    }

    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_cal_he[i][j][0]=0.;
            dssd_cal_he[i][j][1]=1.;
        }
    }

    fcorrcut = -1;
    fmultcut = 10000;
    fisgzstream = false;
    fisranking = true;
    aidaunpkg = new AIDAUnpacker;

    flocalaidaBETA = new AIDA;
    flocalaidaION = new AIDA;
    flocalaidaCORR = new AIDA;
}

BuildAIDAEvents::~BuildAIDAEvents()
{
    //! error on this
    //delete aidaunpkg;

}

void BuildAIDAEvents::Init(char* aidafile)
{
    flastevent = -1;

    //! new stuff
    fflag_firsthit = true;
    flastts = 0;
    fflag_ision =false;

    fflag_addFirstIonHit=false;
    fflag_addFirstBetaHit=false;
    fflag_addFirstHit =false;

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

    //if (fthresholdfile!=NULL) aidaunpkg->read_threshold_table(fthresholdfile);
    read_threshold_table(fthresholdfile);

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

void BuildAIDAEvents::ReadHECalibTable()
{

    //clean up
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_cal_he[i][j][0]=0.;
            dssd_cal_he[i][j][1]=1.;
        }
    }

    ifstream inpf(fcalibfile_he);
    if (inpf.fail()){
        cout<<"No Calibration table for high energy is given"<<endl;
        return;
    }

    cout<<"Start reading calibration table for high energy: "<<fcalibfile_he<<endl;
    Int_t dssd_index,strip_index;
    Double_t cal1,cal2;
    Int_t mm=0;

    while (inpf.good()){
    //for (Int_t i=0;i<100;i++){
        inpf>>dssd_index>>strip_index>>cal1>>cal2;
        dssd_cal_he[dssd_index][strip_index][0]=cal1;
        dssd_cal_he[dssd_index][strip_index][1]=cal2;
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
    //hit->SetEnergy((double)aidaraw.adcData*dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][0]);
    hit->SetEnergy((double)aidaraw.adcData*dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][0]);

    hit->SetTimestamp(aidaraw.extTimestamp*tm_stp_scaler_ratio);
    //hit->SetTimestamp(aidaraw.timestamp);
    hit->SetID(aidaraw.stripNo + aidaraw.dssdNo*NumStrXY);
    hit->SetXY(aidaraw.stripNo);
    hit->SetZ(aidaraw.dssdNo);
    hit->SetFEE(aidaraw.feeNo);
    hit->SetFEEChannel(aidaraw.chNo);
    hit->SetRange(aidaraw.rangeType);
    flocalaidaION->AddHit(hit);
}

//! Add AIDA Beta hitsflocalbeta
void BuildAIDAEvents::AddAIDABetaHits(rawaida_info aidaraw){
    if (aidaraw.adcData>dssd_thr[aidaraw.dssdNo][aidaraw.stripNo]){
        AIDAHit* hit = new AIDAHit;
        hit->SetADC(aidaraw.adcData);
        hit->SetEnergy((double)aidaraw.adcData*dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][0]);
        hit->SetTimestamp(aidaraw.extTimestamp*tm_stp_scaler_ratio);
        //hit->SetTimestamp(aidaraw.timestamp);
        if ((flastfastts[aidaraw.feeNo][aidaraw.chNo] + fwindowDisc) > aidaraw.timestamp)
            hit->SetFastTimestamp(flastfasttsEXT[aidaraw.feeNo][aidaraw.chNo]);
        hit->SetID(aidaraw.stripNo + aidaraw.dssdNo*NumStrXY);
        hit->SetXY(aidaraw.stripNo);
        hit->SetZ(aidaraw.dssdNo);
        hit->SetFEE(aidaraw.feeNo);
        hit->SetFEEChannel(aidaraw.chNo);
        hit->SetRange(aidaraw.rangeType);
        flocalaidaBETA->AddHit(hit);
    }
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


    if (mult<64&&!fflag_trans&&(hitx[0]+hity[0]<fmultcut)&&(hitx[1]+hity[1]<fmultcut)&&(hitx[2]+hity[2]<fmultcut)&&(hitx[3]+hity[3]<fmultcut)&&(hitx[4]+hity[4]<fmultcut)&&(hitx[5]+hity[5]<fmultcut))
    //if (mult<64&&!fflag_trans)
    {
        if (fisranking){//! newly added
            if (flocalaidaBETA->BetaGetPosNew(fcorrcut,fsumexcut,fsumeycut)){
                flocalaidaBETA->SetTimestamp(flocalaidaBETA->GetHit(0)->GetTimestamp());
                if (fflag_filldata) fmtrBETA->Fill();
                fADBetaEntry++;
                return true;
            }else {
                return false;
            }
        }else{
            if (flocalaidaBETA->BetaGetPosAllNew(fcorrcut,fsumexcut,fsumeycut)){
                flocalaidaBETA->SetTimestamp(flocalaidaBETA->GetHit(0)->GetTimestamp());
                if (fflag_filldata) fmtrBETA->Fill();
                fADBetaEntry++;
                return true;
            }else {
                return false;
            }
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

    if (fflag_addFirstHit) {
        flocalaidaION->Clear();
        flocalaidaBETA->Clear();
        if (aidaraw.rangeType==1){
            flastADIonts = aidaraw.timestamp;
            AddAIDAIonHits(aidaraw);
            fflag_ision = true;
        }else if (aidaraw.rangeType==0){
            flastADBetats = aidaraw.timestamp;
            AddAIDABetaHits(aidaraw);
            fflag_ision = false;
        }
    }
    fflag_addFirstHit=false;

    while (!flag_stop){
        if (!aidaunpkg->GetNextHit()) {
            flocalaidaION->Clear();
            flocalaidaBETA->Clear();
            return false;
        }
        aidaraw= aidaunpkg->GetAIDAraw();
        fADcurblock = aidaunpkg->GetCurrentBlock();

        if (aidaraw.infoCode==0){
            if (fflag_firsthit){
                flocalaidaION->Clear();
                flocalaidaBETA->Clear();
                if (aidaraw.rangeType==1){
                    flastADIonts = aidaraw.timestamp;

                }
                else if (aidaraw.rangeType==0){
                    flastADBetats = aidaraw.timestamp;
                }
                flastts = aidaraw.timestamp;
                fflag_firsthit = false;
            }
            if (aidaraw.timestamp<flastts+fwindowHits){
                //! if ion hit
                if (aidaraw.rangeType==1){
                    flastADIonts = aidaraw.timestamp;
                    AddAIDAIonHits(aidaraw);
                    fflag_ision = true;
                }
                //! if beta hit
                else if (aidaraw.rangeType==0){
                    flastADBetats = aidaraw.timestamp;
                    AddAIDABetaHits(aidaraw);
                }
            }else{ //out of window, next event
                if (fflag_ision) {
                    flag_stop = CloseIonEvent();
                    fADtsION = aidaraw.extTimestamp*tm_stp_scaler_ratio;
                    if (!flag_stop) {
                        flocalaidaION->Clear();
                        flocalaidaBETA->Clear();
                        if (aidaraw.rangeType==1){
                            flastADIonts = aidaraw.timestamp;
                            AddAIDAIonHits(aidaraw);
                            fflag_ision = true;
                        }else if (aidaraw.rangeType==0){
                            flastADBetats = aidaraw.timestamp;
                            AddAIDABetaHits(aidaraw);
                            fflag_ision = false;
                        }
                    }else{
                        fflag_addFirstHit=true;
                        fisbeta = false;
                    }
                }else{
                    flag_stop = CloseBetaEvent();
                    fADtsBETA = aidaraw.extTimestamp*tm_stp_scaler_ratio;
                    if (!flag_stop) {
                        flocalaidaION->Clear();//!caution later one can remove it
                        flocalaidaBETA->Clear();
                        if (aidaraw.rangeType==1){
                            flastADIonts = aidaraw.timestamp;
                            AddAIDAIonHits(aidaraw);
                            fflag_ision = true;
                        }else if (aidaraw.rangeType==0){
                            flastADBetats = aidaraw.timestamp;
                            AddAIDABetaHits(aidaraw);
                            fflag_ision = false;
                        }
                    }else{
                        fflag_addFirstHit=true;
                        fisbeta = true;
                    }
                }
            }
            flastts = aidaraw.timestamp;
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

        //! handle correlation scaler
        if (fflag_corrscaler_in_stream&&aidaraw.infoCode==9){
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

    }// flag stop
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

void BuildAIDAEvents::read_threshold_table(char* inf)
{
    //clean up
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_thr[i][j]=-10000.;
        }
    }
    ifstream inpf(inf);
    if (inpf.fail()){
        cout<<"No Threshold file given!"<<endl;
        return;
    }
    cout<<"Start reading threshold table from "<<inf<<endl;
    Int_t dssd_index,strip_index;
    Double_t threshold;
    Int_t mm=0;
    while (inpf.good()){
    //for (Int_t i=0;i<100;i++){
        inpf>>dssd_index>>strip_index>>threshold;
        dssd_thr[dssd_index][strip_index]=threshold;
        //cout<<dssd_thr[dssd_index][strip_index]<<endl;
        mm++;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    flag_threhold=true;
    inpf.close();
}



