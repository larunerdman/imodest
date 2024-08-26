#!/bin/bash

#PBS -l walltime=1:00:00,vmem=2G
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-Batch
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/OE-Batch

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/
JOB=PCA-script-stage1-1000g
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR

CANCERS=(BLCA CESC ESCA KIRC LIHC LUSC PAAD PRAD SARC STAD THCA UCEC BRCA COAD GBM HNSC KIRP LGG LUAD OV PCPG READ SKCM TGCT THYM)

for i in "${CANCERS[@]}" ; do 

	CANCER=$i
	DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/
	#DATA_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/$CANCER/genotype-calls-july-2016/

	OE=$DATA_DIR/OE-$JOB/
	mkdir -p $OE
	rm $OE/*

	qsub $WORKING_DIR/$SH_SCRIPT -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=$CANCER-postqc,OE=$OE \
		 -l walltime=23:59:59,vmem=10G \
		 -o $OE -e $OE \
		 
done 