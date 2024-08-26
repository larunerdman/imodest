#!/bin/bash

#PBS -l vmem=20G,walltime=23:59:59,nodes=1:ppn=6

module load shapeit/2.790

(	>&2 echo chr: $chr
	>&2 echo INPUT: $INPUT
	>&2 echo DATA_DIR: $DATA_DIR
	>&2 echo WORKING_DIR: $WORKING_DIR 
	)

cd $DATA_DIR

shapeit --input-bed $INPUT-chr$chr.bed $INPUT-chr$chr.bim $INPUT-chr$chr.fam \
		--input-map /hpf/projects/arnold/Reference-files-phase3/genetic_map_chr${chr}_combined_b37.txt \
		--output-max $INPUT-chr$chr.haps $INPUT-chr$chr.sample \
		--thread 6 --states 200 

### IMPUTE BY chrOMOSOME
JOB=impute-outerloop-v2
SH_SCRIPT=$JOB.sh

OE=$DATA_DIR/OE-$JOB/chr-$chr/
mkdir -p $OE
rm $OE/*

qsub -v DATA_DIR=$DATA_DIR,INPUT=$INPUT,chr=$chr,WORKING_DIR=$WORKING_DIR \
	 -o $OE -e $OE \
	 -l vmem=2G,walltime=72:00:00,nodes=1:ppn=1,gres=localhd:0 \
	 $WORKING_DIR/$SH_SCRIPT

