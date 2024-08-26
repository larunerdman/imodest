#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00
umask 00007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/clinical_preprocess/
INPUT_DIR=$HOME_DIR/clinical_preprocess_input/
OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/clinical_preprocess_output/

JOB=10_clinical_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R

#### Cancer types

CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)


OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi


for i in {0..28}; do
	# Input data
	PICKED_CLINICAL_DATA=$HOME_DIR/clinical_preprocess_input/${CANCERS[i]}.clin.picked.txt
	#SNP_PCA=$HOME_DIR/SNP_preprocess_input/${CANCERS[i]}/genotype-calls-july-2016/${CANCERS[i]}-postqc-ur-ols-rmvd-pruned-PCS.txt
	#SNP_PCA=$INPUT_DIR/all-TCGA-postqc-ur-ols-rmvd-pruned-PCS.txt
	SNP_PCA=$INPUT_DIR/all-TCGA-mj-sub-postqc-ur-ols-rmvd-pruned-PCS.txt

	#### Output
	mkdir -p $OUTPUT_DIR
	CLINICAL_DATA=$OUTPUT_DIR/${CANCERS[i]}_clinical_data.rds
	GENERAL_CLINICAL_MATRIX=$OUTPUT_DIR/${CANCERS[i]}_general_clinical_matrix.rds
	TUMOR_CLINICAL_MATRIX=$OUTPUT_DIR/${CANCERS[i]}_tumor_clinical_matrix.rds
	NORMAL_CLINICAL_MATRIX=$OUTPUT_DIR/${CANCERS[i]}_normal_clinical_matrix.rds

	R_JOB="module load R/3.3.2;

		   Rscript --vanilla $R_SCRIPT
		   --CANCER_TYPE ${CANCERS[i]}
		   --PICKED_CLINICAL_DATA $PICKED_CLINICAL_DATA
		   --SNP_PCA $SNP_PCA
		   --CLINICAL_DATA $CLINICAL_DATA
		   --GENERAL_CLINICAL_MATRIX $GENERAL_CLINICAL_MATRIX
		   --TUMOR_CLINICAL_MATRIX $TUMOR_CLINICAL_MATRIX
		   --NORMAL_CLINICAL_MATRIX $NORMAL_CLINICAL_MATRIX"

	echo $R_JOB | qsub -N ${JOB} \
					   -l vmem=10G,walltime=1:00:00 \
					   -o $OE -e $OE
done
