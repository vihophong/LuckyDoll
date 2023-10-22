#include "BuildAIDAEventsNew.h"
BuildAIDAEvents::BuildAIDAEvents()
{
    fwindowHits = 201;
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
        fsumexcut[i] = -10000.;
        fsumeycut[i] = -10000.;
        fsumexcuth[i] = -10000.;
        fsumeycuth[i] = -10000.;
    }
    for (Int_t i=0;i<NumDSSD;i++){
        for (Int_t j=0;j<NumStrXY;j++){
            dssd_cal_he[i][j][0]=0.;
            dssd_cal_he[i][j][1]=1.;
        }
    }

    for (int i=0;i<NumFee;i++) for (int j=0;j<NumChFee;j++) chMask[i][j]=1;


    fcorrcut = -1;
    fmultcut = 10000;
    fisgzstream = false;
    fisranking = true;
    aidaunpkg = new AIDAUnpacker;

    flocalaidaBETA = new AIDA;
    flocalaidaION = new AIDA;
    flocalaidaCORR = new AIDA;
    ftprevcheck=-1000;
    //h1=new TH1F("tt","tt",2000,0,2000);
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

    aidaunpkg->CopyChannelMask(chMask);
    aidaunpkg->ResetChannelMask();

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
    //! Add on 2017 May20th

    //! Add on 2017 July 26
    memset(fmultxyz,0,sizeof(fmultxyz));

    //! Added July29, 2017
    memset(foverflowflag,0,sizeof(foverflowflag));

    //! stuff for asics timestamp correction

    fflagasicstscorr=true;
    for (int i=0;i<200;i++){
        fprev_ASICS_cnt[i]=0;
        ffirst_ASICS_ts[i]=0;
        ffirst_ASICS_extts[i]=0;
        fprev_ASICS_ts[i]=0;
        fprev_ASICS_extts[i]=0;
    }

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
    //! stuff for asics timestamp correction
    int asicsNo=aidaraw.feeNo*4+aidaraw.chNo/16;
    unsigned long long corrts=aidaraw.timestamp;
    unsigned long long corrextts=aidaraw.extTimestamp;
    if ((aidaraw.timestamp-fprev_ASICS_ts[asicsNo])==200){
        if (fprev_ASICS_cnt[asicsNo]==0) {
            ffirst_ASICS_ts[asicsNo]=fprev_ASICS_ts[asicsNo];
            ffirst_ASICS_extts[asicsNo]=fprev_ASICS_extts[asicsNo];
        }
        if (fflagasicstscorr){
            corrts=ffirst_ASICS_ts[asicsNo];
            corrextts=ffirst_ASICS_extts[asicsNo];
            //corrts=aidaraw.timestamp- (200*fprev_ASICS_cnt[asicsNo]);
            //corrextts=aidaraw.extTimestamp- (200*fprev_ASICS_cnt[asicsNo]);
        }
        fprev_ASICS_cnt[asicsNo]++;
    }else{
        fprev_ASICS_cnt[asicsNo]=0;
    }
    fprev_ASICS_ts[asicsNo]=aidaraw.timestamp;
    fprev_ASICS_extts[asicsNo]=aidaraw.extTimestamp;
    aidaraw.timestamp=corrts;
    aidaraw.extTimestamp=corrextts;

    //! flag check (for briken2015)
    if (aidaraw.dssdNo<0||aidaraw.stripNo<0) return;
    //if (chMask[aidaraw.feeNo][aidaraw.chNo]==1){ //no disable high energy
        AIDAHit* hit = new AIDAHit;
        hit->SetADC(aidaraw.adcData);
        //hit->SetEnergy((double)aidaraw.adcData*dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][0]);
        hit->SetEnergy((double)aidaraw.adcData*dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal_he[aidaraw.dssdNo][aidaraw.stripNo][0]);
        hit->SetTimestamp(aidaraw.extTimestamp);
        //if ((flastfastts[aidaraw.feeNo][aidaraw.chNo] + fwindowDisc) > aidaraw.timestamp)
        //   hit->SetFastTimestamp(flastfasttsEXT[aidaraw.feeNo][aidaraw.chNo]);
        hit->SetFastTimestamp(aidaraw.timestamp);
        hit->SetID(aidaraw.stripNo + aidaraw.dssdNo*NumStrXY);
        hit->SetXY(aidaraw.stripNo);
        hit->SetZ(aidaraw.dssdNo);
        hit->SetFEE(aidaraw.feeNo);
        hit->SetFEEChannel(aidaraw.chNo);
        hit->SetRange(aidaraw.rangeType);
        flocalaidaION->AddHit(hit);
    //}
}

//! Add AIDA Beta hitsflocalbeta
void BuildAIDAEvents::AddAIDABetaHits(rawaida_info aidaraw){
    //! stuff for asics timestamp correction
    int asicsNo=aidaraw.feeNo*4+aidaraw.chNo/16;
    unsigned long long corrts=aidaraw.timestamp;
    unsigned long long corrextts=aidaraw.extTimestamp;
    if ((aidaraw.timestamp-fprev_ASICS_ts[asicsNo])==200){
        if (fprev_ASICS_cnt[asicsNo]==0) {
            ffirst_ASICS_ts[asicsNo]=fprev_ASICS_ts[asicsNo];
            ffirst_ASICS_extts[asicsNo]=fprev_ASICS_extts[asicsNo];
        }
        if (fflagasicstscorr){
            corrts=ffirst_ASICS_ts[asicsNo];
            corrextts=ffirst_ASICS_extts[asicsNo];
            //corrts=aidaraw.timestamp- (200*fprev_ASICS_cnt[asicsNo]);
            //corrextts=aidaraw.extTimestamp- (200*fprev_ASICS_cnt[asicsNo]);
        }
        fprev_ASICS_cnt[asicsNo]++;
    }else{
        fprev_ASICS_cnt[asicsNo]=0;
    }
    fprev_ASICS_ts[asicsNo]=aidaraw.timestamp;
    fprev_ASICS_extts[asicsNo]=aidaraw.extTimestamp;
    aidaraw.timestamp=corrts;
    aidaraw.extTimestamp=corrextts;

    //! flag check (for briken2015)
    if (aidaraw.dssdNo<0||aidaraw.stripNo<0) return;
    //! make flag
    fmultxyz[aidaraw.dssdNo*128+aidaraw.stripNo]++;
    if(fmultxyz[aidaraw.dssdNo*128+aidaraw.stripNo]>1){
        flocalaidaBETA->AddStripMult(aidaraw.dssdNo);
        if (aidaraw.stripNo<64) flocalaidaBETA->SetStripMultX1FlagMask(aidaraw.dssdNo,aidaraw.stripNo);
        else if (aidaraw.stripNo>=64&&aidaraw.stripNo<128) flocalaidaBETA->SetStripMultX2FlagMask(aidaraw.dssdNo,aidaraw.stripNo);
        else if (aidaraw.stripNo>=128&&aidaraw.stripNo<192) flocalaidaBETA->SetStripMultY1FlagMask(aidaraw.dssdNo,aidaraw.stripNo);
        else if (aidaraw.stripNo>=192) flocalaidaBETA->SetStripMultY2FlagMask(aidaraw.dssdNo,aidaraw.stripNo);
    }
    if (aidaraw.adcData==32768) foverflowflag[aidaraw.dssdNo]++;

    if (aidaraw.adcData>dssd_thr[aidaraw.dssdNo][aidaraw.stripNo]&&chMask[aidaraw.feeNo][aidaraw.chNo]==1){
        AIDAHit* hit = new AIDAHit;
        hit->SetADC(aidaraw.adcData);
        hit->SetEnergy((double)aidaraw.adcData*dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][1] + dssd_cal[aidaraw.dssdNo][aidaraw.stripNo][0]);
        hit->SetTimestamp(aidaraw.extTimestamp);
        //if ((flastfastts[aidaraw.feeNo][aidaraw.chNo] + fwindowDisc) > aidaraw.timestamp)
        //   hit->SetFastTimestamp(flastfasttsEXT[aidaraw.feeNo][aidaraw.chNo]);
        hit->SetFastTimestamp(aidaraw.timestamp);
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


    /*
    cout<<"\n****begin event****"<<endl;
    for (int i=0;i<flocalaidaBETA->GetHits().size();i++){
        int asicsNoall=flocalaidaBETA->GetHit(i)->GetFEE()*4+flocalaidaBETA->GetHit(i)->GetFEEChannel()/16;
        int asicsNo=flocalaidaBETA->GetHit(i)->GetFEEChannel()/16;
        if (asicsNo==2&&flocalaidaBETA->GetHit(i)->GetFEE()==9) cout<<fADIonEntry<<"beta: "<<asicsNoall<<"-"<<flocalaidaBETA->GetHit(i)->GetFEEChannel()<<"-"<<asicsNo<<" original ts="<<flocalaidaBETA->GetHit(i)->GetFastTimestamp()<<" | corrected ts= "<<flocalaidaBETA->GetHit(i)->GetTimestamp()<<endl;
    }
    for (int i=0;i<flocalaidaION->GetHits().size();i++){
        int asicsNo=flocalaidaION->GetHit(i)->GetFEEChannel()/16;
        if (asicsNo==2&&flocalaidaION->GetHit(i)->GetFEE()==9)   cout<<fADIonEntry<<"ion: "<<flocalaidaION->GetHit(i)->GetFEEChannel()<<"-"<<asicsNo<<" original ts="<<flocalaidaION->GetHit(i)->GetFastTimestamp()<<" | corrected ts= "<<flocalaidaION->GetHit(i)->GetTimestamp()<<endl;
    }
    cout<<"\n****end event****"<<endl;
    */

    if (flocalaidaION->GetMult()<=0) return false;
    flocalaidaION->SetTimestamp(flocalaidaION->GetHit(0)->GetTimestamp());
    flocalaidaION->SetPrevIonTimestamp(flastADIonts);
    if (flastADBetats>flastADIonts)
        flocalaidaION->SetTimeWindow((Long64_t)flastADBetats-(Long64_t)fADtsION);
    else
        flocalaidaION->SetTimeWindow((Long64_t)flastADIonts-(Long64_t)fADtsION);


    if (fisranking){//! newly added
        //if (flocalaidaION->IonGetPos()) {
        if (flocalaidaION->IonGetPosNew(-1,fsumexcuth,fsumeycuth)) {
            /*
            for (int i=0;i<flocalaidaBETA->GetHits().size();i++){
                AIDAHit* hit=new AIDAHit;
                flocalaidaBETA->GetHit(i)->Copy(*hit);
                flocalaidaION->AddHit(hit);
            }
            */
            //if (fflag_filldata) fmtrION->Fill();
            //fADIonEntry++;
            //return true;
        }

    }else{
        //if (flocalaidaION->IonGetPos()) {
        if (flocalaidaION->IonGetPosAllNew2(-1,fsumexcuth,fsumeycuth)) {
            /*
            for (int i=0;i<flocalaidaBETA->GetHits().size();i++){
                AIDAHit* hit=new AIDAHit;
                flocalaidaBETA->GetHit(i)->Copy(*hit);
                flocalaidaION->AddHit(hit);
            }
            */
            //if (fflag_filldata) fmtrION->Fill();
            //fADIonEntry++;
            //return true;
        }
    }

    //! check overflow channels
    unsigned short dmaxz=0;

    if (flocalaidaION->GetMaxZ()!=NumDSSD-1){
        for (int i=flocalaidaION->GetMaxZ()+1;i<NumDSSD;i++){
           if (foverflowflag[i]>0||flocalaidaION->GetMultX(i)>0||flocalaidaION->GetMultY(i)>0) dmaxz=(unsigned short)(i-flocalaidaION->GetMaxZ());
        }
    }
    flocalaidaION->SetDeltaMaxZ(dmaxz);   

    //if (fflag_filldata) fmtrION->Fill();
    fADIonEntry++;
    return true;
}

//!close event BETA
bool BuildAIDAEvents::CloseBetaEvent()
{
    if(fverbose>0)
      cout << __PRETTY_FUNCTION__ << endl;

    if (flocalaidaBETA->GetMult()<=0) return false;


    unsigned short* hitx= flocalaidaBETA->GetMultXs();
    unsigned short* hity= flocalaidaBETA->GetMultYs();
    unsigned short mult=flocalaidaBETA->GetMult();
    flocalaidaBETA->SetTimestamp(flocalaidaBETA->GetHit(0)->GetTimestamp());

    //! added july 26,2017
    flocalaidaBETA->SetTimeWindow((Long64_t)flastADBetats-(Long64_t)fADtsBETA);
    flocalaidaBETA->SetPrevIonTimestamp(flastADIonts);

    //if (!fflag_trans&&(hitx[0]+hity[0]<fmultcut)&&(hitx[1]+hity[1]<fmultcut)&&(hitx[2]+hity[2]<fmultcut)&&(hitx[3]+hity[3]<fmultcut)&&(hitx[4]+hity[4]<fmultcut)&&(hitx[5]+hity[5]<fmultcut))
    //if (mult<64&&!fflag_trans&&(hitx[0]+hity[0]<fmultcut)&&(hitx[1]+hity[1]<fmultcut)&&(hitx[2]+hity[2]<fmultcut)&&(hitx[3]+hity[3]<fmultcut)&&(hitx[4]+hity[4]<fmultcut)&&(hitx[5]+hity[5]<fmultcut))
    //if (mult<64&&!fflag_trans)
    //{
        if (fisranking){//! newly added
            if (flocalaidaBETA->BetaGetPosNew(fcorrcut,fsumexcut,fsumeycut)){
                //if (fflag_filldata) fmtrBETA->Fill();
                //fADBetaEntry++;
                //return true;
            }else {
                //return false;
            }
        }else{
            //flocalaidaBETA->BetaGetPosYonly();
            if (flocalaidaBETA->BetaGetPosAllNew2(fcorrcut,fsumexcut,fsumeycut)){
                //if (fflag_filldata) fmtrBETA->Fill();
                //fADBetaEntry++;
                //return true;
            }else {
                //return false;
            }
        }
    //}
    if (mult>=64){
        //flocalaidaBETA->SetTimestamp(flocalaidaBETA->GetHit(0)->GetTimestamp());
        //if (fflag_filldata) fmtrPULSER->Fill();
        fADPulserEntry++;
        if (fflag_pulser_in_stream) return true;
        else return false;
        //!will not give to the "nextevent" scheme!
    }
    //if (fflag_filldata) fmtrBETA->Fill();
    fADBetaEntry++;
    return true;
}

//!AIDA event builder also
bool BuildAIDAEvents::GetNextEvent(){
    //! read from aida
    bool flag_stop=false;

    if (fflag_addFirstHit) {
        flocalaidaION->Clear();
        flocalaidaBETA->Clear();
        if (aidaraw.rangeType==1){
            fADtsION = aidaraw.extTimestamp;
            flastADIonts = aidaraw.extTimestamp;
            AddAIDAIonHits(aidaraw);
            fflag_ision = true;
        }else if (aidaraw.rangeType==0){
            fADtsBETA = aidaraw.extTimestamp;
            flastADBetats = aidaraw.extTimestamp;
            AddAIDABetaHits(aidaraw);
            fflag_ision = false;
        }
        //h1->Fill(aidaraw.timestamp-ftprevcheck);
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
                fADtsION = aidaraw.extTimestamp;
                fADtsBETA = aidaraw.extTimestamp;
                if (aidaraw.rangeType==1){
                    flastADIonts = aidaraw.extTimestamp;
                }
                else if (aidaraw.rangeType==0){
                    flastADBetats = aidaraw.extTimestamp;
                }
                flastts = aidaraw.timestamp;
                fflag_firsthit = false;
            }
            if (aidaraw.timestamp<=flastts+fwindowHits){
                //! if ion hit
                if (aidaraw.rangeType==1){
                    flastADIonts = aidaraw.extTimestamp;
                    AddAIDAIonHits(aidaraw);
                    fflag_ision = true;
                }
                //! if beta hit
                else if (aidaraw.rangeType==0){
                    flastADBetats = aidaraw.extTimestamp;
                    AddAIDABetaHits(aidaraw);
                }
                //h1->Fill(aidaraw.timestamp-ftprevcheck);
            }else{ //out of window, next event
                ftprevcheck=-1000;
                if (fflag_ision) {
                    flag_stop = CloseIonEvent();
                    if (!flag_stop) {
                        flocalaidaION->Clear();
                        flocalaidaBETA->Clear();
                        if (aidaraw.rangeType==1){
                            flastADIonts = aidaraw.extTimestamp;
                            AddAIDAIonHits(aidaraw);
                            fflag_ision = true;
                        }else if (aidaraw.rangeType==0){
                            flastADBetats = aidaraw.extTimestamp;
                            AddAIDABetaHits(aidaraw);
                            fflag_ision = false;
                        }
                    }else{
                        fflag_addFirstHit=true;
                        fisbeta = false;
                    }
                }else{
                    flag_stop = CloseBetaEvent();
                    if (!flag_stop) {
                        flocalaidaION->Clear();//!caution later one can remove it
                        flocalaidaBETA->Clear();
                        if (aidaraw.rangeType==1){
                            flastADIonts = aidaraw.extTimestamp;
                            AddAIDAIonHits(aidaraw);
                            fflag_ision = true;
                        }else if (aidaraw.rangeType==0){
                            flastADBetats = aidaraw.extTimestamp;
                            AddAIDABetaHits(aidaraw);
                            fflag_ision = false;
                        }
                    }else{
                        fflag_addFirstHit=true;
                        fisbeta = true;
                    }
                }
                fADtsION = aidaraw.extTimestamp;
                fADtsBETA = aidaraw.extTimestamp;
            }
            //ftprevcheck=aidaraw.timestamp;
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
                        flastfasttsEXT[feen][chn] = aidaraw.extTimestamp;
                        //cout<<std::dec<<feen<<"-"<<chn<<endl;
                    }else{
                        cout<<"Somethings wrong with the fast discriminator data!"<<endl;
                    }
                }
            }
        }

        //! handle correlation scaler
        if (fflag_corrscaler_in_stream&&aidaraw.rangeType==-1){
            AIDAHit* hit = new AIDAHit;
            //if we dont set all things -> memory leak
            hit->SetADC(0);
            hit->SetEnergy(0.);
            hit->SetTimestamp(aidaraw.extTimestamp);
            hit->SetID(9999);
            hit->SetXY(0);
            hit->SetZ(0);
            hit->SetFEE(0);
            hit->SetFEEChannel(0);
            hit->SetRange(0);
            flocalaidaCORR->AddHit(hit);
        }

    }// flag stop


    //! Clear mult
    memset(fmultxyz,0,sizeof(fmultxyz));
    memset(foverflowflag,0,sizeof(foverflowflag));


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
