#!/bin/bash

umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/RNA_preprocess/
EMAIL=mingjie.mai@sickkids.ca

JOB=08_RNA_preprocess
SH_SCRIPT=$WORKING_DIR/$JOB.sh
QRESUB=1

NUM_CORES=8

CANCERS=(LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV)


for i in {18..21}; do


    CANCER=${CANCERS[i]}

    OE=$WORKING_DIR/OE_$JOB/$CANCER/
    mkdir -p $OE
    if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

    # Input specification
    INPUT_DIR=$HOME_DIR/rawdata/gdc/FASTQ/$CANCER/

    # Output specification
    OUTPUT_DIR=$HOME_DIR/RNA_preprocess_input/$CANCER/exon/
    mkdir -p $OUTPUT_DIR
    #if [ "$(ls -A $OUTPUT_DIR)" ]; then rm $OUTPUT_DIR/*; fi


    for fastq_dir in $INPUT_DIR/*/; do
        qsub -N $CANCER-$JOB \
             -l mem=24G,vmem=24G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:60G \
             -p -512 \
             -o $OE -e $OE \
             -v fastq_dir=$fastq_dir,HOME_DIR=$HOME_DIR,WORKING_DIR=$WORKING_DIR,OUTPUT_DIR=$OUTPUT_DIR,NUM_CORES=$NUM_CORES,QRESUB=$QRESUB \
             $SH_SCRIPT

        qsub -N $CANCER-$JOB \
             -l mem=24G,vmem=24G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:60G \
             -p -512 \
             -o $OE -e $OE \
             -v fastq_dir=$fastq_dir,HOME_DIR=$HOME_DIR,WORKING_DIR=$WORKING_DIR,OUTPUT_DIR=$OUTPUT_DIR,NUM_CORES=$NUM_CORES,QRESUB=$QRESUB \
             $SH_SCRIPT
    done
done






