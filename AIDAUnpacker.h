#ifndef AIDAUNPACKER_H
#define AIDAUNPACKER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "TH1.h"

#include "rawaida.h"
#include "AIDAdefs.h"

#define HEADER_SIZE 0x18 // Size of header in bytes
#define BLOCK_SIZE 0x10000 //Max block size is 64kb. Amount of useful data given in header
#define MAIN_SIZE 0xFFE8 //main_size = block_size - header_size

using namespace std;

class AIDAUnpacker
{
public:
    AIDAUnpacker();
    virtual ~AIDAUnpacker();

    //! Initialize the input file
    void Init(char* inputfile);

    //! Get first block (must be call)
    void ReadHeader();

    //!Set verbose level
    void SetVerbose(int verbose){fverbose = verbose;}

    //! Set current block
    void SetCurrentBlock(int currentblock){fcurrentblk=currentblock;}

    //! Set current package
    void SetCurrentPackage(int currentpackage){fcurrentpkg=currentpackage;}

    //! List of available FEE card
    bool enableFEEs[33];

    //! Get current block number
    int GetCurrentBlock(){return fcurrentblk;}

    //! Get number of block
    int GetNBlock(){return fnblock;}

    //! Get aida raw data
    rawaida_info GetAIDAraw(){return rawaida;}

    //! Get midas raw data
    rawmidas_info GetMIDASraw(){return midas;}

    //! Get the first syncronization pulse
    int GetFirstSync();

    //! Event loop : return full aida data
    bool GetNextHit();

    //! Get number of good hits
    Long64_t GetHitNumber(){return rawaida.evt;}

    //! Read mapping
    void read_mapping(char* mapping_file);

    //! Read threhold table
    void read_threshold_table(char* inf);

    //! Get correlation scaler histogram
    TH1F* GetCorrScalerHisto(){return corrTS;}

    //! Book a tree
    int BookTree();
    //! Get a tree
    TTree* GetTree(){return rawtree;}

private:
    //! file path
    char* ffilepath;
    //! file name
    ifstream finfile;
    //! verbose level
    int fverbose;
    //! Correlation scaler histogram
    TH1F* corrTS;
    TTree* rawtree;

    //! First Scan flag
    bool first_scan_flag;

    //! number of blocks
    int fnblock;
    //! indicate current block of data
    int fcurrentblk;
    //! indicate current package of data in a block
    int fcurrentpkg;
    //! indicate current data length in a block
    unsigned int fcurrentlen;
    //! end of block flag
    bool feobflag;

    //! container of block data
    char fblkData[MAIN_SIZE];
    char fblkHeader[HEADER_SIZE];

    //! raw midas container
    rawmidas_info midas;

    //! raw aida container
    rawaida_info rawaida;


    //! some unnessary info
    unsigned short int header_stream;
    unsigned short int header_tape;
    unsigned short int header_MyEndian;
    unsigned short int header_DataEndian;

    //! first syncronization flag and correlation scaler
    int first_sync_flag[NumFee];
    int first_corr_scaler_datum;

    //! first correlation scaler
    long long first_corr_scaler_timestamp;
    //! first correlation scaler offset
    long long my_first_time_offset;
    //! first sync on all enabled modules
    unsigned long first_tm_stp_msb_modules[NumFee];

    //!global time stamp MSB
    unsigned long my_tm_stp_msb;
    //!local time stamp MSB
    unsigned long tm_stp_msb_modules[NumFee];
    //! global sync flag
    bool sync_flag;
    //! local sync flag
    bool sync_flag_modules[NumFee]; //true= last SYNC100 received correctly
    //! global pause flag
    bool pause_flag;
    //! local pause flag
    bool pause_flag_modules[NumFee]; //false= synchronization is running (i.e. not paused)

    //! flag for mbs(external timestamp data) hit in each fee
    bool MBS_hit[NumFee][2];
    //! LSB bit of MBS data
    long long MBS_tm_stp_lsb[NumFee][2];
    //! contain a piece of external time stamp data
    long long MBS_bits[NumFee][2];
    long long my_time_offset;


    //! global previous time stamp
    long long tm_stp_prev;
    //! local previous time stamp
    long long my_tm_stp_prev[NumFee]; //check only for adc...? or all

    //! prev external time stamp
    bool fail_ext_tmstmp_Flag;
    long long my_ext_timestamp_prev;

    //!Mapping
    bool flag_mapping;
    int FEEtoDSSD[NumFee][NumChFee];
    int FEEtoStrip[NumFee][NumChFee];
    int chMask[NumFee][NumChFee]; //newly added
    long long feeMask;

    //! Threshold table
    bool flag_threhold;
    Double_t dssd_thr[NumDSSD][NumStrXY];

    //! fill flag for sorter
    bool fillFlag;

    //! Event loop : return the midas data
    bool GetNextHitMIDAS();

    //! Reconstruct raw aida data
    bool ReconstructRawAIDA();

    //! Helper for MBS time stamp lookup
    bool checkCoincidence(const unsigned long t1, const unsigned long t2, const unsigned long t3, unsigned long dt);

    void ClearSorter();
};

#endif // AIDAUNPACKER_H
