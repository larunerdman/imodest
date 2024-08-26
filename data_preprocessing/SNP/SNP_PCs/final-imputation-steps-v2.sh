#!/bin/bash


(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo INPUT: $INPUT
>&2 echo DATA_DIR: $DATA_DIR)


cd $DATA_DIR

mkdir -p imputation-info-files

rm ./imputation-info-files/*
mv $INPUT-chr*-int*_* ./imputation-info-files/

### CREATING GEN FILES BY CATTING INT FILES

#chr=15
#for chr in {12,16,18,19} ; do
for chr in {1..22} ; do 

	cat $INPUT-chr$chr-int* > $INPUT-chr$chr.gen ; 

done 

#chr=15
#for chr in {12,16,18,19} ; do
for chr in {1..22} ; do 

	sed 's/-9/0/' $INPUT-chr$chr.sample > $INPUT-revised-chr$chr.sample ; 
	
done

JOB=gen-qc-v3
SH_SCRIPT=$JOB.sh

### POST-IMPUTATION QC BY chrOMOSOME
module purge
#chr=15
#for chr in {12,16,18,19} ; do
OE=$DATA_DIR/OE-$JOB/
mkdir -p $OE
rm $OE/*

for chr in {1..22} ; do 

	qsub -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=$INPUT,chr=$chr \
		 -l walltime=23:59:59,vmem=40G \
		 -o $OE -e $OE \
		 $WORKING_DIR/$SH_SCRIPT

	sleep 0.01
done
