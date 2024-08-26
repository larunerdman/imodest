#!/bin/bash

#PBS -l walltime=1:00:00,vmem=2G
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-Batch
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-Batch

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/
JOB=final-imputation-steps-v2
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR

#CANCERS=(BLCA CESC ESCA KIRC LIHC LUSC PAAD PRAD SARC STAD THCA UCEC BRCA COAD GBM HNSC KIRP LGG LUAD OV PCPG READ SKCM TGCT)
# CANCERS=(COAD LIHC)
# CANCERS=(BLCA CESC ESCA GBM HNSC)
CANCERS=(KIRC LUSC PAAD PRAD SARC STAD THCA UCEC BRCA KIRP LGG LUAD OV PCPG READ SKCM TGCT BLCA CESC ESCA GBM HNSC)

for i in "${CANCERS[@]}" ; do 
	CANCER=$i
	DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/imputation
	#DATA_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/$CANCER/genotype-calls-july-2016/imputation/

	OE=$DATA_DIR/OE-$JOB/
	mkdir -p $OE
	rm $OE/*


	qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=$CANCER-postqc2 \
		 -l walltime=23:59:59,vmem=8G \
		 -o $OE -e $OE \
		 $SH_SCRIPT
done 