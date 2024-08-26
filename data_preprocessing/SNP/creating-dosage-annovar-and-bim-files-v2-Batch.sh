#!/bin/bash

umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/SNP_preprocess/

JOB=creating-dosage-annovar-and-bim-files-v2
SH_SCRIPT=$WORKING_DIR/$JOB.sh
CANCERS=(LIHC BLCA KIRC PAAD SARC THCA BRCA ESCA KIRP LUAD PCPG SKCM THYM CESC GBM LAML LUSC PRAD STAD UCEC COAD HNSC LGG OV READ TGCT)

for i in {1..25}; do

    CANCER=${CANCERS[i]}

    OE=$WORKING_DIR/OE_$JOB/$CANCER/
    if [ -d "$OE)" ]; then rm -r $OE; fi
    mkdir -p $OE


    DATA_DIR=$HOME_DIR/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/imputation/

    SAMPLE_FILE=$DATA_DIR/$CANCER.sample
    cat $DATA_DIR/$CANCER-postqc2-revised-chr1.sample | grep "TCGA" | cut -c 1-12 > $SAMPLE_FILE

    for j in {1..22} ; do

        diff $DATA_DIR/$CANCER-postqc2-revised-chr$j.sample $DATA_DIR/$CANCER-postqc2-revised-chr1.sample

        qsub -N $CANCER-$JOB-$j \
             -v chr=$j,WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,CANCER=$CANCER \
             -l vmem=20G,walltime=23:59:59,gres=localhd:0 \
             -p 512 \
             -o $OE -e $OE \
             $SH_SCRIPT
    done

done
