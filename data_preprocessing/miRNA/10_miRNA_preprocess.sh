#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/miRNA_preprocess/

JOB=10_miRNA_preprocess
R_SCRIPT=$WORKING_DIR/${JOB}.R

OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

## Reference data
MIRNA_TARGET=$HOME_DIR/ref/targetScan.v7.1/Summary_Counts.all_transcripts.txt
MIRNA_FAMILY=$HOME_DIR/ref/targetScan.v7.1/miR_Family_Info.txt
GENE_GTF=$HOME_DIR/ref/gencode.v19/gencode.v19.annotation.gtf

## Output
mkdir -p $HOME_DIR/ref/target_matrix/
MIRNA_TARGET_MATRIX_GENE=$HOME_DIR/ref/target_matrix/10_miRNA_target_matrix_gene.rds
MIRNA_TARGET_MATRIX_EXON=$HOME_DIR/ref/target_matrix/10_miRNA_target_matrix_exon.rds

R_JOB="module load R/3.3.2;

       Rscript --vanilla $R_SCRIPT
               --MIRNA_TARGET $MIRNA_TARGET
               --MIRNA_FAMILY $MIRNA_FAMILY
               --GENE_GTF $GENE_GTF
               --MIRNA_TARGET_MATRIX_GENE $MIRNA_TARGET_MATRIX_GENE
               --MIRNA_TARGET_MATRIX_EXON $MIRNA_TARGET_MATRIX_EXON"

echo $R_JOB | qsub -N ${JOB} \
                   -l vmem=60G,walltime=23:59:59 \
                   -o $OE -e $OE
