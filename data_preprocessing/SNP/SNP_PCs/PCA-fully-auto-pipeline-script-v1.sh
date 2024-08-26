#!/bin/bash

module load R 

# “d” and “datadir” have the required argument of the directory the data is located in
# “i” and “input” have the required argument of the input file

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

## Get king robust PCs with 1000g reference set
$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT -r 

## Use PCAiR to ID outliers relative to 1000g reference set
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-hap-ceu-tsi-merged-pruned -r

## Remove 1000g-based outliers and get king robust PCs of data alone
$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT-ur -e $INPUT-ur-hapmap-outliers-to-rmv.txt

## Use PCAiR to ID outliers relative to our given set
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-ur-pruned

## Remove outliers relative to our given set and run king robust on remaining individuals
$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT-ur-ols-rmvd -e $INPUT-ur-outliers-to-rmv.txt

## Calculate final PCs
Rscript --vanilla $WORKING_DIR/PCA-pcair-script-v1.R -d $DATA_DIR -i $INPUT-ur-ols-rmvd-pruned -p
