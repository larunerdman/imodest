#!/bin/bash

#PBS -l mem=20G,vmem=20G,walltime=23:00:00

working_dir=/hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/COAD-READ/

setwd $working_dir

module load qctool

awk '{if(length($4) == 1 && length($5) == 1) print}' $working_dir/READ-COAD-postqc2-chr$chr.gen > $working_dir/READ-COAD-no2ab-chr$chr.gen

qctool -g $working_dir/READ-COAD-no2ab-chr$chr.gen \
	   -og $working_dir/READ-COAD-postqc2-postqc-chr$chr.gen \
	   -snp-missing-rate 0.05 \
	   -maf 0.01 1 \
	   -info 0.8 1 \
	   -hwe 20


# /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/converting-imputed-results-to-plink-step1.sh $DATA_DIR $INPUT $chr
