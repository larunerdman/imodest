#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/lncRNA_preprocess/

JOB=40_lncRNA_preprocess
OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
rm $OE/*

TEMP_DIR=$WORKING_DIR/temp/
mkdir -p $TEMP_DIR
rm $TEMP_DIR/*

R_SCRIPT=$WORKING_DIR/$JOB.R


#### Constant, Reference files
NUM_NUCLEOTIDES=1000000
UCSC_EXON_INFO=$HOME_DIR/cnv_preprocess/Ucsc_Exon_Info.rds
LNCRNA_DATA=$WORKING_DIR/gencode.v19.long_noncoding_RNAs.gtf
OUTPUT_DIR=$WORKING_DIR

EXON_SPECIFIC=0

R_JOB="module load R/3.3.2; \
       Rscript --vanilla $R_SCRIPT 
       --NUM_NUCLEOTIDES $NUM_NUCLEOTIDES \
       --EXON_SPECIFIC $EXON_SPECIFIC \
       --UCSC_EXON_INFO $UCSC_EXON_INFO \
       --LNCRNA_DATA $LNCRNA_DATA \
       --OUTPUT_DIR $OUTPUT_DIR \
       --TEMP_DIR $TEMP_DIR"

echo $R_JOB | qsub -N ${JOB}-${EXON_SPECIFIC} \
                       -l vmem=30G,walltime=23:59:59 \
                       -o $OE -e $OE

for j in {1..24}; do
    
    c=$j

    EXON_SPECIFIC=1

    R_JOB="module load R/3.3.2; \
           Rscript --vanilla $R_SCRIPT 
           --NUM_NUCLEOTIDES $NUM_NUCLEOTIDES \
           --EXON_SPECIFIC $EXON_SPECIFIC \
           --UCSC_EXON_INFO $UCSC_EXON_INFO \
           --LNCRNA_DATA $LNCRNA_DATA \
           --OUTPUT_DIR $OUTPUT_DIR \
           --c $c \
           --TEMP_DIR $TEMP_DIR"
   
    echo $R_JOB | qsub -N ${JOB}-${EXON_SPECIFIC}-$c \
                       -l vmem=30G,walltime=23:59:59 \
                       -o $OE -e $OE
done
