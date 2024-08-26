#!/bin/bash

CANCERS=(BLCA LAML PAAD SKCM BRCA LGG PCPG CESC ESCA LIHC PRAD STAD GBM LUAD TGCT HNSC LUSC READ-COAD THCA KIRC SARC THYM KIRP OV UCEC)
#CANCERS=BLCA
#CANCERS=(PAAD SKCM BRCA LGG PCPG CESC ESCA LIHC PRAD STAD GBM LUAD TGCT HNSC LUSC READ-COAD THCA KIRC SARC THYM KIRP OV UCEC)

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/
 
cd $WORKING_DIR
 
for i in "${CANCERS[@]}" ; do 
	CANCER=$i
	DATA_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/data/individual_cohort/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/imputation
	INPUT=$CANCER-cancer-postqc-postqc2
		
	qsub final-imputation-steps-v2.sh -v WORKING_DIR=$WORKING_DIR,DATA_DIR=$DATA_DIR,INPUT=$INPUT ; sleep 3 ; 
	
done 
