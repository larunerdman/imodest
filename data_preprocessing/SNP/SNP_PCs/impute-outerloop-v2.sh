#!/bin/bash


(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo DATA_DIR: $DATA_DIR
>&2 echo INPUT: $INPUT
>&2 echo chr: $chr)

cd $DATA_DIR

IFS=$'\n'

JOB=impute-innerloop-v2
SH_SCRIPT=$JOB.sh

OE=$DATA_DIR/OE-$JOB/chr-$chr/
mkdir -p $OE
rm $OE/*

#for chr in {12,16,18,19} ; do
#for chr in $( seq 1 22 ) ; do
#chr=15	
while read line ; do
	
	NUM_JOB=`qstat -u lauren | awk '{print $10}' | grep "[Q|R]" | wc -l`
	
	while [[ $NUM_JOB -gt 5000 ]]; do
		sleep 1m
		NUM_JOB=`qstat -u lauren | awk '{print $10}' | grep "[Q|R]" | wc -l`
	done

	qsub -v LINE=$line,chr=$chr,INPUT=$INPUT,DATA_DIR=$DATA_DIR,WORKING_DIR=$WORKING_DIR \
 		 -o $OE -e $OE \
		 -l vmem=40G,walltime=23:59:00,gres=localhd:0 \
		 $WORKING_DIR/$SH_SCRIPT

	sleep 0.1

done < int-chr$chr.txt
