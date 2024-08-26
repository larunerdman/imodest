#!/bin/bash

#PBS -l mem=20G,vmem=20G,walltime=23:00:00


(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo DATA_DIR: $DATA_DIR
>&2 echo INPUT: $INPUT
>&2 echo chr: $chr)


cd $DATA_DIR

module load gtool

awk '{if(length($5) == 1 && length($6) == 1) print}' $INPUT-chr$chr.gen > $INPUT-no2ab-chr$chr.gen

grep -vwF -f $INPUT-chr$chr-somaticvars.txt $INPUT-no2ab-chr$chr.gen > $INPUT-no2ab-nosom-chr$chr.gen

qctool -S --g $INPUT-no2ab-nosom-chr$chr.gen \
	   --og $INPUT-postqc-chr$chr.gen \
	   --missing 0.05 \
	   -threshold 0.8


# /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/converting-imputed-results-to-plink-step1.sh $DATA_DIR $INPUT $chr
