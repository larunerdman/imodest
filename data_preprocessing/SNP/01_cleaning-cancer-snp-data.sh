#!/bin/bash

## PREPARING THE CANCER SNP DATA FOR IMPUTATION

DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/data/individual_cohort/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/
WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo CANCER: $CANCER

module load plink/1.90b3x
module load R 

cd /hpf/largeprojects/agoldenb/mingjie/projects/MoR/data/individual_cohort/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/



if [ "$CANCER" == "READ-COAD" ] ; then 

	Rscript --vanilla /hpf/largeprojects/agoldenb/lauren/Modes/scripts/data-processing/get-cancer-inds.R --CANCER $CANCER	
	plink --bfile $CANCER-no- --keep $CANCER-cancer-inds.txt --make-bed --out $CANCER-cancer --noweb
	
else 
	
	cp $CANCER-no-.fam $CANCER-no-.tfam
	Rscript --vanilla /hpf/largeprojects/agoldenb/lauren/Modes/scripts/data-processing/get-cancer-inds.R --CANCER $CANCER
	
	plink --tfile $CANCER-no- --keep $CANCER-cancer-inds.txt --make-bed --out $CANCER-cancer --noweb

fi 

JOB='03-calling-and-preprocessing-TCGA-SNPs'
SH_SCRIPT=$JOB.sh
OE=$WORKING_DIR/OE-$JOB/
mkdir -p $OE

$WORKING_DIR/preprocessing-cancer-snp-data-Aug082018.sh $CANCER-cancer

plink --bfile $CANCER-cancer-postqc --make-bed --out $CANCER-cancer-postqc-postqc

qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,CANCER=$CANCER-cancer-postqc -l walltime=48:00:00,mem=20g,vmem=20g -o $OE -e $OE $WORKING_DIR/$SH_SCRIPT





