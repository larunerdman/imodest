#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/SNP_preprocess/
JOB=making-tped-file

OE=$WORKING_DIR/OE-$JOB/
mkdir -p $OE

CANCERS=(LIHC BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP OV PRAD SKCM THCA)

for i in {1..25}; do

CANCER=${CANCERS[i]}

INPUT_DIR=$HOME_DIR/SNP_preprocess_input/${CANCERS[i]}/genotype-calls-july-2016/

SH_JOB="module load python/2.7.9; \

python $WORKING_DIR/making-tped-file6.py \
-c $INPUT_DIR/birdseed-v2.calls.txt \
-a $WORKING_DIR/GenomeWideSNP_6.na35.annot.csv \
-o $INPUT_DIR/$CANCER.tped"

echo $SH_JOB | qsub -N $CANCER-$JOB \
-l walltime=150:00:00,vmem=40G \
-o $OE -e $OE
done 