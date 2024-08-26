#!/bin/bash

#PBS -l walltime=23:00:00,vmem=2G
#PBS -o /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/OE-Batch
#PBS -e /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/OE-Batch

WORKING_DIR=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess/geno-calling-files/imputation/
JOB=imputation-script-take2
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR


CANCERS=(BRCA THCA KIRC LGG LUAD)


for i in {0..0}; do

	CANCER=${CANCERS[i]}
	DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/$CANCER/genotype-calls-july-2016/imputation/
	mkdir -p $DATA_DIR
	rm $DATA_DIR/$CANCER-postqc2*
	cp $DATA_DIR/../$CANCER-postqc2* $DATA_DIR
	#DATA_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/$CANCER/genotype-calls-july-2016/imputation/

	OE=$DATA_DIR/OE-$JOB/
	mkdir -p $OE
	rm $OE/*

	qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=$CANCER-postqc2 \
		 -l walltime=23:00:00,vmem=10G \
		 -o $OE -e $OE \
		 $SH_SCRIPT
done 