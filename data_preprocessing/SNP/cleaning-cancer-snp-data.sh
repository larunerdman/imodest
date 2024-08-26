#!/bin/bash

module load plink/1.90b3x
module load R 

/hpf/largeprojects/agoldenb/mingjie/projects/MoR/data/individual_cohort/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/

#Rscript --vanilla /hpf/largeprojects/agoldenb/lauren/Modes/scripts/data-processing/get-cancer-inds.R --CANCER $CANCER

plink --tfile $CANCER-no- --remove $CANCER-postqc-healthy-inds.txt --make-bed --out $CANCER-cancer-postqc --noweb

DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/data/individual_cohort/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/
WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files
JOB='03-calling-and-preprocessing-TCGA-SNPs'
SH_SCRIPT=$JOB.sh
OE=$WORKING_DIR/OE-$JOB/
mkdir -p $OE

qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,CANCER=all-TCGA -l walltime=48:00:00,mem=20g,vmem=20g -o $OE -e $OE $SH_SCRIPT





