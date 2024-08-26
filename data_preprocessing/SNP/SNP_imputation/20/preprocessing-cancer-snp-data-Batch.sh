#!/bin/bash

#PBS -l walltime=1:00:00,vmem=2G
#PBS -o /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/OE-Batch
#PBS -e /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/OE-Batch

WORKING_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/
JOB=preprocessing-cancer-snp-data
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR


OE=$WORKING_DIR/OE-$JOB/
mkdir -p $OE
rm $OE/*


CANCERS=(BRCA THCA KIRC LGG LUAD)


for i in {0..0}; do

	CANCER=${CANCERS[i]}
	DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/$CANCER/genotype-calls-july-2016/
	#DATA_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/$CANCER/genotype-calls-july-2016/

	qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,CANCER=$CANCER \
		 -l walltime=23:59:59,vmem=30G \
		 -o $OE -e $OE \
		 $SH_SCRIPT
done 