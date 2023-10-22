#!/bin/bash                                                                

# make file list                                                  
#ls -tr /home/aida/raid/aida/rawdata/takechi_exp/R$1* | wc -l > list/list_R$1.txt
#ls -tr /home/aida/raid/aida/rawdata/takechi_exp/R$1* >> list/list_R$1.txt

if [ $# == 4 ]; then

    expr $3 - $2 + 1 > list/list_R$4_$1_$2to$3.txt
    for i in `seq $2 $3`;
    do
	ls -tr aidaraw/R$1_$i >> list/list_R$4_$1_$2to$3.txt
    done  
    # run
    ./aidapartial -a list/list_R$4_$1_$2to$3.txt -o rootfiles/$4_aida$1_$2to$3.root -map config/FEE_table_briken2016.txt -thr config/thr_r141_1cps.txt -cal config/cal_briken2016_r25_26.txt -wi 5000

else
    echo "usage: ./aidabuild_takechi.sh run_number subrun_number_low subrun_number_high yymmdd_hhmm"
fi




