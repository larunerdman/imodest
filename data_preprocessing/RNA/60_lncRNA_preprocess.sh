#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/lncRNA_preprocess/

JOB=60_lncRNA_preprocess
R_SCRIPT=$WORKING_DIR/${JOB}.R

OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

## Constants
NUM_NUCLEOTIDES=1000000

## Reference data
GENE_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf
LNCRNA_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.long_noncoding_RNAs.gtf

## Output
TARGET_MATRIX_DIR=$HOME_DIR/ref/target_matrix/
mkdir -p $TARGET_MATRIX_DIR
LNCRNA_TARGET_MATRIX_GENE=$TARGET_MATRIX_DIR/60_lncRNA_target_matrix_gene.rds
LNCRNA_TARGET_MATRIX_EXON=$TARGET_MATRIX_DIR/60_lncRNA_target_matrix_exon.rds

R_JOB="module load R/3.3.2;

       Rscript --vanilla $R_SCRIPT
               --NUM_NUCLEOTIDES $NUM_NUCLEOTIDES
               --GENE_GTF $GENE_GTF
               --LNCRNA_GTF $LNCRNA_GTF
               --LNCRNA_TARGET_MATRIX_GENE $LNCRNA_TARGET_MATRIX_GENE
               --LNCRNA_TARGET_MATRIX_EXON $LNCRNA_TARGET_MATRIX_EXON"

echo $R_JOB | qsub -N ${JOB} \
                   -l vmem=60G,walltime=48:00:00 \
                   -o $OE -e $OE
