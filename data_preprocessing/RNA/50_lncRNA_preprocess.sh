#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=${HOME_DIR}/lncRNA_preprocess/

JOB=50_lncRNA_preprocess
OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi


R_SCRIPT=${WORKING_DIR}/${JOB}.R

INPUT_DIR=${WORKING_DIR}/temp/
OUTPUT_DIR=${HOME_DIR}/lncRNA_preprocess_output/

for j in {1..1}; do 
    EXON_SPECIFIC=$j

    R_JOB="module load R/3.3.2; \
       Rscript --vanilla $R_SCRIPT 
       --INPUT_DIR $INPUT_DIR \
       --EXON_SPECIFIC $EXON_SPECIFIC \
       --OUTPUT_DIR $OUTPUT_DIR"

    echo $R_JOB | qsub -N ${JOB} \
                       -l vmem=30G,walltime=23:59:59 \
                       -o $OE -e $OE
done

