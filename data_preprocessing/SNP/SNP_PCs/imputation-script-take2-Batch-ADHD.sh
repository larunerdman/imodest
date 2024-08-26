#!/bin/bash

#PBS -l walltime=23:00:00,vmem=2G
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-Batch
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-Batch

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/
JOB=imputation-script-take2
SH_SCRIPT=$JOB.sh

cd $WORKING_DIR


DATAFILE=(psychchip hcore snp6)

for i in {0..2}; do

#	CANCER=${CANCERS[i]}
	DATA_DIR=/hpf/largeprojects/agoldenb/lauren/schachar-collab-2017/${DATAFILE[i]}/imputation/
	mkdir -p $DATA_DIR
#	rm $DATA_DIR/$CANCER-postqc2*
#	cp $DATA_DIR/../$CANCER-postqc2* $DATA_DIR
	#DATA_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/$CANCER/genotype-calls-july-2016/imputation/

	OE=$DATA_DIR/OE-$JOB/
	mkdir -p $OE
	rm $OE/*

	qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=${DATAFILE[i]}-postqc \
		 -l walltime=23:00:00,vmem=10G \
		 -o $OE -e $OE \
		 $SH_SCRIPT
done 