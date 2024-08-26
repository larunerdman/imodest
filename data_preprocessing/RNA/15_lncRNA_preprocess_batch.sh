#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/lncRNA_preprocess/

JOB=15_lncRNA_preprocess
R_SCRIPT=$WORKING_DIR/${JOB}.R
CANCERS=(BLCA COAD HNSC LAML LUAD PAAD READ STAD THYM BRCA ESCA KIRC LGG LUSC PCPG SARC TGCT UCEC CESC GBM KIRP LIHC OV PRAD SKCM THCA)

OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

for i in {0..25}; do
	CANCER_TYPE=${CANCERS[i]}
	RNA_EXPRESSION_MATRIX=$HOME_DIR/RNA_preprocess_output/${CANCERS[i]}_RNASeq_gene.rds
	LNCRNA_IDS=$WORKING_DIR/10_lncRNA_ids.rds
	LNCRNA_EXPRESSION_MATRIX=$HOME_DIR/lncRNA_preprocess_output/${CANCERS[i]}_lncRNA_gene_expression.rds

	R_JOB="module load R/3.3.2; \
		   Rscript --vanilla $R_SCRIPT 
		   --CANCER_TYPE ${CANCERS[i]} \
		   --RNA_EXPRESSION_MATRIX $RNA_EXPRESSION_MATRIX \
		   --LNCRNA_IDS $LNCRNA_IDS \
		   --LNCRNA_EXPRESSION_MATRIX $LNCRNA_EXPRESSION_MATRIX"
	if [ -f $LNCRNA_EXPRESSION_PROFILE ]; then
		echo $R_JOB | qsub -N $CANCER_TYPE-${JOB} \
						   -l vmem=10G,walltime=2:00:00 \
						   -o $OE -e $OE
	fi
	
done

