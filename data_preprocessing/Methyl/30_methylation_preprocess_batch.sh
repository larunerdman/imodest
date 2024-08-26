#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/methylation_preprocess/

JOB=30_methylation_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R
OE=$WORKING_DIR/OE_$JOB/

# Input data
CANCERS=(READ-COAD BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP LIHC OV PRAD SKCM THCA)

mkdir -p $OE
rm $OE/*

for i in {0..0}; do

	CANCER_TYPE=${CANCERS[i]}
	METHYLATION_INPUT_FOLDER=$HOME_DIR/methylation_preprocess_input/
	METHYLATION_OUTPUT_FOLDER=$HOME_DIR/methylation_preprocess_output/raw/

		R_JOB="	module load R/3.3.2; \
				Rscript --vanilla $R_SCRIPT \
				--CANCER_TYPE $CANCER_TYPE \
				--METHYLATION_INPUT_FOLDER $METHYLATION_INPUT_FOLDER \
				--METHYLATION_OUTPUT_FOLDER $METHYLATION_OUTPUT_FOLDER"

		echo $R_JOB | qsub -N $CANCER_TYPE-$JOB \
						   -l vmem=60G,walltime=2:00:00 \
						   -e $OE -o $OE
done
