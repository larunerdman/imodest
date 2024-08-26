#!/bin/bash

#PBS -l walltime=23:00:00,vmem=2G
#PBS -o /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-PCA-pipeline-script-v1
#PBS -e /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/imputation/OE-PCA-pipeline-script-v1

# “d” and “datadir” have the required argument of the directory the data is located in
# “i” and “input” have the required argument of the input file

module load R 
module load plink

# read the options
TEMP=`getopt -o d:i: --long datadir:,input: -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -d|--datadir)
            case "$2" in
                *) DATA_DIR=$2 ; shift 2 ;;
            esac ;; 
		-i|--input)
            case "$2" in
                *) INPUT=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files
hapmap_dir=/hpf/largeprojects/agoldenb/lauren/Reference-data-1000g-PCA/

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

cd $DATA_DIR


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
echo Creating final \$CANCER-postqc set
echo
echo =================================================
plink --bfile $INPUT-ur-ols-rmvd-ols-rmvd --make-bed --out ${INPUT}2 