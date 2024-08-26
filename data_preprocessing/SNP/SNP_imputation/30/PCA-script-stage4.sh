#!/bin/bash

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

cd $DATA_DIR

module load eigensoft


twstats -t /hpf/tools/centos6/eigensoft/6.0/src/smarttables/twtable \
		-i ${INPUT}2-pca.eigenvalue \
		-o ${INPUT}2-pca.twstats