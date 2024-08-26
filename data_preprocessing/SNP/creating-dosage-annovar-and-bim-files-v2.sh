#!/bin/bash

umask 0007

module load python/2.7.9

INPUT=$DATA_DIR/$CANCER-postqc2-postqc-chr${chr}

python $WORKING_DIR/recoding-gen-files.py -g $INPUT.gen -o $INPUT.dosage

awk -v chr=$chr '{print chr,$3,$4,0,$5,$6}' $INPUT.gen > $INPUT.bim


rm $DATA_DIR/$CANCER-postqc2-postqc-chr${chr}.gen
