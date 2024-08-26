#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/methylation_preprocess/

JOB=40_methylation_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R

OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

TEMP_DIR=$WORKING_DIR/temp/
mkdir -p $TEMP_DIR
if [ "$(ls -A $TEMP_DIR)" ]; then rm $TEMP_DIR/*; fi

NUM_NUCLEOTIDES=1000000
GENE_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf
METHYLATION_PROBE_INFO_450=$WORKING_DIR/probeInfo_450.rds


for j in {1..24}; do

    c=$j

    for k in {0..0}; do

        EXON_SPECIFIC=$k

        R_JOB=" module load R/3.3.2;

                Rscript --vanilla $R_SCRIPT
                --NUM_NUCLEOTIDES $NUM_NUCLEOTIDES
                --c $c
                --EXON_SPECIFIC $EXON_SPECIFIC
                --GENE_GTF $GENE_GTF
                --METHYLATION_PROBE_INFO_450 $METHYLATION_PROBE_INFO_450
                --TEMP_DIR $TEMP_DIR"

        echo $R_JOB |   qsub -N ${JOB}-${EXON_SPECIFIC} \
                             -l vmem=40G,walltime=23:59:59 \
                             -o $OE -e $OE
    done
done

