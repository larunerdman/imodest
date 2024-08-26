#!/bin/bash

#PBS -l walltime=23:00:00,vmem=2G
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-PCA-pipeline-script-v1
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-PCA-pipeline-script-v1

module load R 
module load plink

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

cd $DATA_DIR

echo =================================================
echo 
echo Subsetting to samples overlapping MJ\'s set
echo
echo =================================================
Rscript --vanilla $WORKING_DIR/subsetting-to-mj-samples.R -d $DATA_DIR -i $INPUT

plink --bfile $INPUT --keep $INPUT-ol-with-ther-data-keep-ids.txt -make-bed --out $INPUT-ovl-dat

INPUT=$INPUT-ovl-dat

echo =================================================
echo 
echo Getting king robust PCs with 1000g reference set
echo
echo =================================================
$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT -r -e ""

echo =================================================
echo
echo Using PCAiR to ID outliers relative to 1000g reference set
echo
echo =================================================
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-hap-ceu-tsi-merged-pruned -r

echo =================================================
echo
echo Removing 1000g-based outliers and get king robust PCs of data alone
echo
echo =================================================

$WORKING_DIR/PCA-king-script-v1.sh -e $INPUT-hap-ceu-tsi-merged-pruned-hapmap-outliers-to-rmv.txt -d $DATA_DIR -i $INPUT-ur 

echo =================================================
echo
echo Using PCAiR to ID outliers relative to our given set
echo
echo =================================================
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-ur-pruned

echo =================================================
echo
echo Removing outliers relative to our given set and run king robust on remaining individuals
echo
echo =================================================
$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT-ur-ols-rmvd -e $INPUT-ur-pruned-outliers-to-rmv.txt

echo =================================================
echo
echo Calculating final PCs
echo
echo =================================================
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-ur-ols-rmvd-pruned -p

echo =================================================
echo
echo Creating final \$CANCER-postqc2 set
echo
echo =================================================
plink --bfile $INPUT-ur-ols-rmvd-ols-rmvd --make-bed --out ${INPUT}2 

echo =================================================
echo
echo Copying \$CANCER-postqc2.{bed,bim,fam} to ./imputation/ directory
echo
echo =================================================
cp ${INPUT}2.{bed,bim,fam} $DATA_DIR/imputation/