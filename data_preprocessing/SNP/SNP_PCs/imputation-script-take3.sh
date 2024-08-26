#!/bin/bash

(
>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo DATA_DIR: $DATA_DIR
>&2 echo INPUT: $INPUT
)

cd $DATA_DIR

module load plink/1.90b3x
module load R/3.2.3

### REMOVE AMBIGUOUS SNPS

Rscript --vanilla $WORKING_DIR/identifying-ambig-snps.R -i $INPUT
Rscript --vanilla $WORKING_DIR/removing-dup-positions.R -b $INPUT

cat ambig-snps.txt duplicated-position-snps.txt > ambig-and-dup-snps.txt

plink --bfile $INPUT --exclude ambig-and-dup-snps.txt --make-bed --out $INPUT-noambig --noweb

### SPLIT PLINK FILE BY chrOMOSOMES

for i in {1..22} ; do 

	plink --bfile $INPUT-noambig --chr $i --make-bed --out $INPUT-chr$i --noweb ; 

done

### GENERATING INT FILES
for chr in {1..22} ; do

	Rscript --vanilla $WORKING_DIR/int-file-creation.R -i $INPUT -c $chr ; 

done

### PREPHASE BY chrOMOSOME
JOB=shapeit-innerloop-v3
SH_SCRIPT=$JOB.sh
OE=$DATA_DIR/OE-$JOB/
mkdir -p $OE
rm $OE/*

module purge
for chr in {1..22} ; do 

	qsub -v chr=$chr,INPUT=$INPUT,DATA_DIR=$DATA_DIR,WORKING_DIR=$WORKING_DIR \
		 -l vmem=40G,walltime=23:59:00,nodes=1:ppn=8,gres=localhd:0 \
		 -o $OE -e $OE \
		 $WORKING_DIR/$SH_SCRIPT

		 sleep 0.1

done

