#include "AIDAUnpackerGz.h"

AIDAUnpacker::AIDAUnpacker():midas(),rawaida(),finfile(),fifsinfilegz(),finfilegz()
{
    fenableFastDisc = false;
    rawtree=NULL;
    //! Default threshold
    for (int i=0;i<NumDSSD;i++){
        for (int j=0;j<NumStrXY;j++){
            dssd_thr[i][j]=0.;
        }
    }
    fmaxtsoffset = 100;
    fisgzstream = false;
    for (int i=0;i<NumFee;i++) for (int j=0;j<NumChFee;j++) chMask[i][j]=1;
    fncorrscaler=0;
}

AIDAUnpacker::~AIDAUnpacker(){
    //! error on this
    //delete corrTS;
    if (!fisgzstream) finfile.close();
    else fifsinfilegz.close();
}

void AIDAUnpacker::Init(char *inputfile){
    fverbose = 0;
    nresetwarning = 0;
    first_scan_flag=true;
    fillFlag = true;
    //!Set default Fee mask (enable all)
    feeMask = FEEMASK;

    //! Default mapping and threshold table flag
    flag_mapping = false;
    flag_threhold = false;

    //! ReSet first sync pulse for each modules
    for(int i=0;i<NumFee;i++){
        first_sync_flag[i] = 0;
        first_tm_stp_msb_modules[i] = 0;
    }
    //! ReSet first correlation time stamp
    first_corr_scaler_datum = 0;
    first_corr_scaler_timestamp = 0;
    my_first_time_offset = 0;


    if (!fisgzstream){//Normal mode
        //! Open files
        ffilepath = inputfile;
        finfile.open(ffilepath,std::ios::in|std::ios::binary);
        if (!finfile.good()){
            cout<<"Error opening file "<<ffilepath<<endl;
            exit(0);
        }
        //! Check for data integrity and get number of blocks
        finfile.seekg(0,finfile.end);
        int fileSize = finfile.tellg();
        finfile.seekg(0,finfile.beg);
        fileSize += -finfile.tellg();

        fnblock=  fileSize / BLOCK_SIZE;
        int nBlock_check = fileSize % BLOCK_SIZE;
        if (nBlock_check != 0){
            cout<<"Missing block data!"<<endl;
            exit(0);
        }
        cout<<__PRETTY_FUNCTION__<< " Openning raw aida file: "<<ffilepath<<endl;
        cout<< "Number of block = "<<fnblock<<endl;

    }else{//Gzip mode
        ffilepath = inputfile;
        fifsinfilegz.open(ffilepath, std::ios::in|std::ios::binary);
        if (!fifsinfilegz.good()){
            std::cout<<"Error opening file "<<ffilepath<<std::endl;
            exit(0);
        }
        finfilegz.push(boost::iostreams::gzip_decompressor());
        finfilegz.push(fifsinfilegz);
        fnblock= 2000000000;//just put a finite large number
        cout<<__PRETTY_FUNCTION__<< " Openning GZipped raw aida file: "<<ffilepath<<endl;
    }

    //!Read header of the first block
    fcurrentblk = 0;
    ReadHeader();
    if (fisgzstream&&finfilegz.eof()) exit(0);
    ClearSorter();
}

void AIDAUnpacker::ClearSorter(){
    rawaida.evt=0;
    //!Clear SOME LOCAL Var for sorter
    for (int i=0;i<NumFee;i++){
        tm_stp_msb_modules[i]=0;
        sync_flag_modules[i]=false;
        pause_flag_modules[i]=false;
        my_tm_stp_prev[i]=0;
        MBS_hit[i][0]=false;
        MBS_hit[i][1]=false;
        MBS_tm_stp_lsb[i][0]=0;
        MBS_tm_stp_lsb[i][1]=0;
        MBS_bits[i][0]=0;
        MBS_bits[i][1]=0;
        my_tm_stp_msb=0;
    }
    my_time_offset=0;
    tm_stp_prev=0;

}

void AIDAUnpacker::ReadHeader(){
    fcurrentpkg = 0;
    if (!fisgzstream){// read normally
        finfile.read((char*)&fblkHeader,sizeof(fblkHeader));
        finfile.read((char*)&fblkData,sizeof(fblkData));
    }else{//read from gz pipe
        finfilegz.read((char*)&fblkHeader,sizeof(fblkHeader));
        finfilegz.read((char*)&fblkData,sizeof(fblkData));
    }
    //!READ HEADER
    //! some unnessecaray stuffs
    header_stream = (fblkHeader[12] & 0xFF) << 8 | (fblkHeader[13]& 0xFF);
    header_tape = (fblkHeader[14] & 0xFF) << 8 | (fblkHeader[15]& 0xFF);
    header_MyEndian = (fblkHeader[16] & 0xFF) << 8 | (fblkHeader[17]& 0xFF);
    header_DataEndian = (fblkHeader[18] & 0xFF) << 8 | (fblkHeader[19]& 0xFF);
    //! current data length
    fcurrentlen =
    (fblkHeader[20] & 0xFF) | (fblkHeader[21]& 0xFF) << 8 |
    (fblkHeader[22] & 0xFF) << 16  | (fblkHeader[23]& 0xFF) << 24 ;
}

bool AIDAUnpacker::GetNextHitMIDAS(){
    midas.Clear();
    unsigned int word0; //each word contains 32 bit
    unsigned int word1;
    word0 = (fblkData[fcurrentpkg] & 0xFF)  | (fblkData[fcurrentpkg+1] & 0xFF) << 8 |
      (fblkData[fcurrentpkg+2] & 0xFF) << 16 | (fblkData[fcurrentpkg+3] & 0xFF) << 24;
    word1 = (fblkData[fcurrentpkg+4] & 0xFF) | (fblkData[fcurrentpkg+5] & 0xFF) << 8 |
      (fblkData[fcurrentpkg+6] & 0xFF) << 16 | (fblkData[fcurrentpkg+7] & 0xFF) << 24;

    //!Read next block if there is a end of block word : AIDA format: only one FFFFFFFF word in end of block
    if( (word0 & 0xFFFFFFFF) == 0xFFFFFFFF || (word1 & 0xFFFFFFFF) == 0xFFFFFFFF){
        ReadHeader();
        word0 = (fblkData[fcurrentpkg] & 0xFF)  | (fblkData[fcurrentpkg+1] & 0xFF) << 8 |
          (fblkData[fcurrentpkg+2] & 0xFF) << 16 | (fblkData[fcurrentpkg+3] & 0xFF) << 24;
        word1 = (fblkData[fcurrentpkg+4] & 0xFF) | (fblkData[fcurrentpkg+5] & 0xFF) << 8 |
          (fblkData[fcurrentpkg+6] & 0xFF) << 16 | (fblkData[fcurrentpkg+7] & 0xFF) << 24;
        fcurrentblk++;
    }

    //! After getting two word, try to extract data
    midas.dataType = (word0 >> 30) & 0x3; //mask at bit 31+30

    //!Adc data
    if(midas.dataType == 0x3){//ADC data
        unsigned int chIdent = ( word0>>16 ) & 0x00000FFF; //channel Identification on bits 16:27
        midas.feeId = (chIdent>>6) & 0x0000003F;
        midas.chId = chIdent & 0x0000003F;
        midas.adcRange = (word0 >>28) & 0x00000001; //bit 28
        midas.adcData = word0&0x0000FFFF; //bit 0:15
        midas.timestampLsb = word1 & 0x0FFFFFFF; //bit 0:27
    }

    //! Information data needed for getting MSB time stamp,
    //! sync100 syncronization and timestamp colleration scaler
    else if(midas.dataType == 0x2){//information data
        midas.feeId = (word0 >> 24) & 0x0000003F; //bits 24:29
        midas.infoField = word0 & 0x000FFFFF; //bits 0:19
        midas.infoCode = (word0 >> 20) & 0x0000000F; //bits 20:23
        midas.timestampLsb = word1 & 0x0FFFFFFF;  //bits 0:27
    }
    else if(midas.dataType == 0x1){//sample length data
        unsigned int chIdent =  (word0 >> 16) & 0x00000FFF; //bits 16:27
        midas.feeId = (chIdent >> 6) & 0x0000003F;
        midas.chId = chIdent & 0x0000003F;
        midas.sampleLengh = word0 & 0x0000FFFF; //bits 0:15
        midas.timestampLsb = word1 & 0x0FFFFFFF;  //bits 0:27
    }
    else if(midas.dataType == 0x0){//Wave form data (usually follow by sample length data
        midas.waveformtrace1= (word0 >> 16) & 0x00003FFF;
        midas.waveformtrace2= 0x00003FFF;
        midas.waveformtrace3= (word1 >> 16) & 0x00003FFF;
        midas.waveformtrace4= 0x00003FFF;
    }

    if (fverbose > 0) cout << __PRETTY_FUNCTION__ << "current block= "<<fcurrentblk << " / "<< fnblock << "current package= "<< fcurrentpkg << endl;
    if (fverbose > 1) midas.Show();

    fcurrentpkg += 8; // For every 8 byte (64 bit) for a block data
    if (fcurrentpkg>=fcurrentlen){
        ReadHeader();
        //! to fix timejump bug!
        if(!fisgzstream&&finfile.eof()) return false;

        if (fisgzstream&&finfilegz.eof()) return false;
        fcurrentblk++;
    }
    if (fcurrentblk>fnblock) return false;
    else return true;
}

bool AIDAUnpacker::GetNextHit(){
    if (GetNextHitMIDAS()){
        //! reconstruct raw aida data  here
        if (!first_scan_flag){
            while(!ReconstructRawAIDA()){
                if (!GetNextHitMIDAS())
                    return false;
            }
        }else{
            ReconstructRawAIDA();
        }
        if (rawtree!=NULL&&!first_scan_flag) {
            if (ts_prev==rawaida.timestamp&&prev_rt==rawaida.rangeType&&prev_fee==rawaida.feeNo&&prev_ch==rawaida.chNo)
                cout<<"duplicate entry found!"<<rawaida.feeNo <<
                      " cr blk"<<fcurrentblk<<" prev blk"<<prev_blk<<
                      " cr pkg"<<fcurrentpkg<<" prev pkg"<<prev_pkg<<endl;
            ts_prev=rawaida.timestamp;
            prev_rt=rawaida.rangeType;
            prev_fee=rawaida.feeNo;
            prev_ch=rawaida.chNo;
            prev_blk=fcurrentblk;
            prev_pkg=fcurrentpkg;
            rawtree->Fill();
            ncheck++;
        }
        rawaida.evt++;
        return true;
    }else{
        return false;
    }
    //up to here we have valid midas entry
}

bool AIDAUnpacker::checkCoincidence(const unsigned long t1, const unsigned long t2, const unsigned long t3, unsigned long dt)
{
    if( t2 > (t1+dt) ){
      return false;
      //std::cout << "one"<<std::endl;
    }
    if( t1 > (t2+dt) ){
      //std::cout << "two"<<std::endl;
      return false;}
    if( t1 > (t3+dt) ){
      //std::cout << "3one"<<std::endl;
      return false;}
    if( t3 > (t1+dt) ){
      //std::cout << "4one"<<std::endl;
      return false;}

    return true;
}

bool AIDAUnpacker::ReconstructRawAIDA(){
    rawaida.Clear();
    // get status of flags from previous entry
    sync_flag=sync_flag_modules[midas.feeId];
    pause_flag=pause_flag_modules[midas.feeId];

    fillFlag=true;
    if ((feeMask>>midas.feeId)&0x1==0){
        fillFlag=false;
    }
    //AIDA ADC data (just fill in some things,
    //timestamp data for LSB will be used later when there is a fillFlag
    else if (midas.dataType==3){
        rawaida.infoCode=0;
        rawaida.rangeType=midas.adcRange;
        rawaida.feeNo=midas.feeId;
        rawaida.chNo=midas.chId;
        rawaida.adcData=midas.adcData;
    }
    else if (midas.dataType==2&&midas.infoCode!=0){ //GREAT data info ->To get timestamp of MSB, also omit Unidentify data
        //First get aida fee ID for this info
        rawaida.feeNo=midas.feeId;
        rawaida.infoCode=midas.infoCode;

        //!Get first sync
        if (midas.infoCode==2||midas.infoCode==3||midas.infoCode==4){ //get msb of module
            if (first_sync_flag[midas.feeId]==0){ //get first dantum with time stamp msb info
                first_tm_stp_msb_modules[midas.feeId]=(midas.infoField & 0x000FFFFF);
                //std::cout<<"Got you first sync pulse on Fee module No"<<midas.feeId<<" MSB of timestamp= "<<std::hex<<"0x"<<first_tm_stp_msb_modules[midas.feeId]<<std::dec<<std::endl;
            }
            first_sync_flag[midas.feeId]++;
        }

        //if SYNC100 pulse:
        //get tm stp MSB +
        if (midas.infoCode==4){
            my_tm_stp_msb=(midas.infoField & 0x000FFFFF);
            tm_stp_msb_modules[midas.feeId]=my_tm_stp_msb;
            sync_flag_modules[midas.feeId]=true;
            sync_flag=true;
        }
        //if PAUSE time stamp
        //get tm stp MSB +
        else if (midas.infoCode==2){
            my_tm_stp_msb=(midas.infoField & 0x000FFFFF);
            tm_stp_msb_modules[midas.feeId]=my_tm_stp_msb;
            pause_flag_modules[midas.feeId]=true;
            pause_flag=true;
        }
        //if RESUME time stamp
        else if (midas.infoCode==3){
            pause_flag_modules[midas.feeId]=false;
            my_tm_stp_msb=(midas.infoField & 0x000FFFFF);
            tm_stp_msb_modules[midas.feeId]=my_tm_stp_msb;
            pause_flag=false;
        }
        // if discriminator data >>we need to KEEP
        else if (midas.infoCode==6){
            rawaida.chNo=(midas.infoField & 0x000FFFFF);
            rawaida.adcData=(midas.infoField&0x000FFFFF); //newversion of midas firmware
        }
        //if MBS information index
        //-> Reconstruct EXT time from correlation scaler
        //-> Recontstruct AIDA time stamp
        else if (midas.infoCode==8){
            //!omit start bit of time stamp
            rawaida.infoCode = -1;

            int my_MBS_index;
            long long my_MBS_bits;
            my_MBS_index=(midas.infoField & 0x000F0000) >>16; //-> Bit 19:16 is information index
            my_MBS_bits=midas.infoField & 0x0000FFFF; //-> Bit 15:0 with EXT scaler data (0 1 2 index) (LUPO?)


            if (my_MBS_index==0 || my_MBS_index==1){ //Get low bits of EXT time stamp
                MBS_hit[midas.feeId][my_MBS_index]=true;
                MBS_bits[midas.feeId][my_MBS_index]=my_MBS_bits; // LSB of the ext data (bit 15:0 or 31:16) with ch number
                MBS_tm_stp_lsb[midas.feeId][my_MBS_index]=midas.timestampLsb; //Get timestamp of this information
            }
            else if(my_MBS_index==2){ //get higher bits of EXT timestamp (assume other low bits info already come
                fncorrscaler++;
                rawaida.rangeType = -1;// newly added (May13 2017)

                if (MBS_hit[midas.feeId][0]&&MBS_hit[midas.feeId][1]){ // when we have enough 2 hits in time stamp scaler
                    unsigned long t1=midas.timestampLsb;
                    unsigned long t2=MBS_tm_stp_lsb[midas.feeId][1];
                    unsigned long t3=MBS_tm_stp_lsb[midas.feeId][0];

                    //Another condition to check if t1 is in coincidence with t2 t3 within 1 us window (so the information come from 1 burse) (100 pulses)
                    if (checkCoincidence(t1,t2,t3,100)){
                        //!Only take last bit of time stamp
                        rawaida.infoCode = 8;
                        //note my_MBS_bits<<32 is now in index 2 which is bit 32-47
                        rawaida.extTimestamp=(my_MBS_bits <<32 | MBS_bits[midas.feeId][1] <<16 | MBS_bits[midas.feeId][0])*tm_stp_scaler_ratio;

                        //cout<<"phong1 - "<<rawaida.extTimestamp<<endl;
                        //!Get first correlation scaler
                        if (first_sync_flag[midas.feeId]>0){
                            if (first_corr_scaler_datum==0) {
                                first_corr_scaler_timestamp=rawaida.extTimestamp;
                                //get offset time from here!!!
                                long long timestamp_temp=(long long)(tm_stp_msb_modules[midas.feeId] << 28 ) | (midas.timestampLsb & 0x0FFFFFFF); //caution this conversion!
                                my_first_time_offset=first_corr_scaler_timestamp-timestamp_temp;
                                //corrTS=new TH1F("corrTS","corrTS",fmaxtsoffset*2,my_first_time_offset-fmaxtsoffset,my_first_time_offset+fmaxtsoffset);
                                //std::cout<<"Got you first offset bw AIDA-LUPO! = "<<std::dec<<my_first_time_offset<<"-- Timestamp LUPO="<<first_corr_scaler_timestamp<<"| Timestamp AIDA="<<timestamp_temp<<std::endl;
                                //CAUTION:IF FOR SOME REASON WE MISS SYNC THEN THIS OFFSET MAKE NON SENSE
                            }
                            first_corr_scaler_datum++;
                        }

                        //addition condition to imply the condition of monotopicialy increasing time stamp
                        if ((rawaida.extTimestamp>my_ext_timestamp_prev)&&(!fail_ext_tmstmp_Flag)) {
                            fillFlag=true;
                        }else if((rawaida.extTimestamp>my_ext_timestamp_prev)&&(fail_ext_tmstmp_Flag)){ //condition to resume
                            fail_ext_tmstmp_Flag=false; //resume
                            fillFlag=true;
                        }
                        else{
                            fail_ext_tmstmp_Flag=true; //pause
                        }

                        if (!fail_ext_tmstmp_Flag) my_ext_timestamp_prev=rawaida.extTimestamp;
                    }else{
                        std::cout<<"MBS info is OUT of Sync FEE#"<<int(midas.feeId)<<", MBS info(info, time-stamp)= ("
                                << my_MBS_bits<<", "<<midas.timestampLsb<<"), ("
                                << MBS_bits[midas.feeId][1]<<", "<<MBS_tm_stp_lsb[midas.feeId][1]<<"), ("
                                << MBS_bits[midas.feeId][0]<<", "<<MBS_tm_stp_lsb[midas.feeId][0]<<")"
                                <<std::endl;
                    }
                    //reset MBS hit and tmp stamp
                    MBS_hit[midas.feeId][0]= false;
                    MBS_hit[midas.feeId][1]= false;
                    MBS_tm_stp_lsb[midas.feeId][0]=0;
                    MBS_tm_stp_lsb[midas.feeId][1]=0;
                    MBS_bits[midas.feeId][0]=0;
                    MBS_bits[midas.feeId][1]=0;
                }
            }

        }// if (midas.infoCode==8)
    }//if great info data
    //Skip uninteresting piece of data!
    if(fillFlag&&!first_scan_flag){
        //get time stamp which include msb
        rawaida.timestamp=(long long)(tm_stp_msb_modules[rawaida.feeNo] << 28 ) | (midas.timestampLsb & 0x0FFFFFFF); //caution this conversion!
        //check reconstructed timestamp monotonically increasing
        // if not... probably missed some synchronization pulse
        if(rawaida.timestamp < my_tm_stp_prev[rawaida.feeNo]){
            //std::cout<<"OUT of Sync nnaida"<<rawaida.feeNo<<" : ";
            //std::cout<<rawaida.timestamp<<std::hex<<"-"<<my_tm_stp_prev[rawaida.feeNo]<<std::hex<<std::endl;
            sync_flag_modules[rawaida.feeNo]=false;
            sync_flag=false;
        }
        my_tm_stp_prev[rawaida.feeNo]=rawaida.timestamp;

        //!mapping, convert ADC and reconstruct ext timestamp
        if (flag_mapping){
            if (rawaida.chNo<NumStrXY&&rawaida.feeNo<NumFee){
                rawaida.dssdNo=FEEtoDSSD[rawaida.feeNo][rawaida.chNo]; //channel index start from 1 since index must start from 1
                rawaida.stripNo=FEEtoStrip[rawaida.feeNo][rawaida.chNo];

                if (rawaida.stripNo<NumStrX) rawaida.adcData=pow(2,15)-rawaida.adcData;//X strip, Negative polarity
                if (rawaida.stripNo>=NumStrX) rawaida.adcData=rawaida.adcData-pow(2,15);//Y strip, Positive polarity
            }else{
                //must be discriminator data
                //print first within first 1000 event
                //if (rawaida.evt<1000) cout<<"Channel No./Fee No. error!: FeeNo.="<<rawaida.feeNo<<"-ChNo.="<<rawaida.chNo<<"-InfoCode="<<rawaida.infoCode<<endl;
                //sync_flag=false;
            }
        }

        if (sync_flag){
            //! Get correlation scaler offset
            if (rawaida.infoCode==8) {
                long long temp=rawaida.extTimestamp-rawaida.timestamp;
                if (temp-my_first_time_offset<fmaxtsoffset&&temp-my_first_time_offset>-fmaxtsoffset) {
                    my_time_offset=temp;

                }else{
                    if (nresetwarning<10) cout<<"AIDA time stamp scaler jumped! maybe timestamp reset is sent"<<endl;
                    nresetwarning ++;
                }
                //corrTS->Fill(temp);
            }else{
                rawaida.extTimestamp=my_time_offset+rawaida.timestamp; //convert to ext timestamp unit
            }

            //!Check global time warps
            if (tm_stp_prev>rawaida.timestamp){
                cout<<"AIDA time stamp has gone down!"<<std::dec<<tm_stp_prev<<"--"<<rawaida.timestamp<<" fee"<<rawaida.feeNo<<endl;
                fillFlag=false;
            }
            tm_stp_prev=rawaida.timestamp;

            //! Masking for slow discriminator data
            if (rawaida.infoCode==0&&rawaida.rangeType==0){
                if(chMask[rawaida.feeNo][rawaida.chNo]){
                    //if (rawaida.infoCode==0) corrTS->Fill(rawaida.extTimestamp-rawaida.timestamp);
                }else{
                    fillFlag=false;
                }
                //! masking for data less than threshold
                if (flag_threhold&&rawaida.adcData<dssd_thr[rawaida.dssdNo][rawaida.stripNo])
                    fillFlag=false;
            }

            //! Masking for high energy
            if (rawaida.infoCode==0&&rawaida.rangeType==1){
                if(chMask[rawaida.feeNo][rawaida.chNo]){

                }else{
                    fillFlag=false;
                }
            }

        }else{ //fail sync flag
            fillFlag=false;
        }
    }else{ //fail fill flag (from MBS hit scan)
        fillFlag=false;
    }
    if (!(rawaida.infoCode==0||(fenableFastDisc&&rawaida.infoCode==6)||rawaida.rangeType==-1)) return false;
    return fillFlag;
}

void AIDAUnpacker::read_mapping(char* mapping_file){
    std::ifstream inpf(mapping_file,std::ios::in);
    if (inpf.fail()){
        cout<<"No Mapping table is given"<<endl;
        return;
    }
    std::cout<<"Start reading mapping from "<<mapping_file<<std::endl;
    //reset
    for (Int_t i=0;i<NumFee;i++){
        for (Int_t j=0;j<NumChFee;j++){
            FEEtoDSSD[i][j]=-1;
            FEEtoStrip[i][j]=-1;
        }
    }
    //do mapping
    Int_t lineread=0;
    Int_t FEE_index,ch_index,DSSD_index,strip_index,ch_mask;
    //basiclly index=fee_index*ch_index!
    while (inpf.good()){
        inpf>>FEE_index>>ch_index>>DSSD_index>>strip_index>>ch_mask;
        //std::cout<<FEE_index<<ch_index<<DSSD_index<<strip_index<<std::endl;
        FEEtoDSSD[FEE_index][ch_index]=DSSD_index;
        FEEtoStrip[FEE_index][ch_index]=strip_index;

        //! reverse mapping!
        DSSDtoFee[DSSD_index][strip_index]=FEE_index;
        DSSDtoCh[DSSD_index][strip_index]=ch_index;

        chMask[FEE_index][ch_index]=ch_mask;
        lineread++;
    }
    std::cout<<"Read "<<lineread<<" line"<<endl;
    flag_mapping=true;
    inpf.close();
}

void AIDAUnpacker::read_threshold_table(char* inf)
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

int AIDAUnpacker::GetFirstSync(){
    int nhits=0;
    int ntotalhits=0;

    int flag[NumFee];
    for (int i=0;i<NumFee;i++){
        if ((feeMask>>i)&0x1) ntotalhits++;
        flag[i]=0;
    }
    cout<<ntotalhits<<endl;

    while (GetNextHit()){
        if (first_tm_stp_msb_modules[rawaida.feeNo]>0){
            if (((feeMask>>rawaida.feeNo)&0x1)&&flag[rawaida.feeNo]==0){
                nhits++;
                flag[rawaida.feeNo]++;
            }
        }
        if (nhits==ntotalhits&&first_corr_scaler_timestamp>0) break;
    }
    ClearSorter();

    for (int i=0;i<33;i++){
        if (first_tm_stp_msb_modules[i]>0){
            cout<<"FEE No. "<<std::dec<<i;
            cout<<" First MSB of Timestamp = "<<"0x"<<std::hex<<first_tm_stp_msb_modules[i]<<endl;
        }
    }
    cout<<"First Corr timestamp = "<<std::dec<<first_corr_scaler_timestamp<<endl;
    cout<<"First Offset timestamp = "<<my_first_time_offset<<endl;

    //!Set starting value of correlation scaler and most significant bit of timestamp
    for (int i=0;i<33;i++){
        tm_stp_msb_modules[i]= first_tm_stp_msb_modules[i];
    }
    my_time_offset = my_first_time_offset;
    first_scan_flag=false;
    if (!fisgzstream){
        //! Rewind
        finfile.clear();
        finfile.seekg(0, ios::beg);
    }else{
        //! Rewind
        fifsinfilegz.clear();
        fifsinfilegz.seekg(0, ios::beg);
        //! reset boost gz pipe
        finfilegz.reset();
        //! setup boost gz pipe again
        finfilegz.push(boost::iostreams::gzip_decompressor());
        finfilegz.push(fifsinfilegz);
    }
    //!Read header of the first block
    fcurrentblk = 0;
    ReadHeader();
    ncheck=0;

    fncorrscaler=0;
}

int AIDAUnpacker::BookTree(TTree *tree)
{
    rawtree=tree;
    rawtree->Branch("evt",&rawaida.evt,"evt/L");
    rawtree->Branch("timestamp",&rawaida.timestamp,"timestamp/L");
    rawtree->Branch("extTimestamp",&rawaida.extTimestamp,"extTimestamp/L");

    rawtree->Branch("feeNo",&rawaida.feeNo,"feeNo/I");
    rawtree->Branch("chNo",&rawaida.chNo,"chNo/I");
    rawtree->Branch("dssdNo",&rawaida.dssdNo,"dssdNo/I");
    rawtree->Branch("stripNo",&rawaida.stripNo,"stripNo/I");

    rawtree->Branch("infoCode",&rawaida.infoCode,"infoCode/S");
    rawtree->Branch("rangeType",&rawaida.rangeType,"rangeType/S");
    rawtree->Branch("adcData",&rawaida.adcData,"adcData/I");
}
