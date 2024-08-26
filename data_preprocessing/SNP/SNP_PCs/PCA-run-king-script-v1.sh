#!/bin/bash

# “i” and “input” have the required argument of the input file
# “r” and “reference” have no arguments, flags use of reference sample in input data

# set an initial value for the reference data flag
REF=0

# read the options
TEMP=`getopt -o d:re::i: --long datadir:,reference,exclude::,input: -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -d|--datadir)
            case "$2" in
                *) DATA_DIR=$2 ; shift 2 ;;
            esac ;; 
        -r|--reference) REF=1 ; shift ;;
		-e|--exclude)
            case "$2" in
				"") EXCLUDE=0 ; shift 2;;
                *) EXCLUDE=$2 ; shift 2 ;;
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

if [ $REF -eq 1 ]
then 
	echo "Using 1000Genomes Caucasian reference samples"
	$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT -r
else
	if [ $EXCLUDE -eq 0 ]
	then 
		$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT
	else
		$WORKING_DIR/PCA-king-script-v1.sh -d $DATA_DIR -i $INPUT -e $EXCLUDE
	fi
fi

