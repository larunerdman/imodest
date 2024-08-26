#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/SNP_preprocess/

# Script description
JOB=creating-dosage-annovar-and-bim-files-v2
SH_SCRIPT=$WORKING_DIR/$JOB.sh

#CANCERS=(LIHC COAD)
CANCERS=(LIHC COAD BLCA GBM KIRP OV PRAD SKCM THCA BRCA HNSC LAML LUAD PAAD READ STAD THYM CESC ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC)

for i in {0..25}; do

    CANCER=${CANCERS[i]}
    OE=$WORKING_DIR/IO_$JOB/$CANCER/
    mkdir -p $OE
    if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

    DATA_DIR=$HOME_DIR/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/imputation/post-qc-gen-sample-files/
    mkdir -p $DATA_DIR
    if [ "$(ls -A $DATA_DIR)" ]; then rm $DATA_DIR/*; fi

    cp $DATA_DIR/../$CANCER-postqc2-postqc* $DATA_DIR

    SAMPLE_FILE=$DATA_DIR/$CANCER.sample
    cat $DATA_DIR/../$CANCER-postqc2-revised-chr1.sample | grep "TCGA" | cut -c 1-12 > $SAMPLE_FILE

    for j in {1..22}; do

        diff $DATA_DIR/../$CANCER-postqc2-revised-chr$j.sample $DATA_DIR/../$CANCER-postqc2-revised-chr1.sample

        qsub -N $CANCER-$JOB-$j\
             -v chr=$j,WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,CANCER=$CANCER \
             -l vmem=20G,walltime=23:59:59,gres=localhd:0 \
             -p 512 \
             -o $OE -e $OE \
             $SH_SCRIPT
    done

done
