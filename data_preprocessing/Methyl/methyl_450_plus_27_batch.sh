#!/bin/bash

#PBS -l vmem=1G,walltime=01:00:00

MY_HOME=/hpf/largeprojects/agoldenb/mingjie
WORKING_DIR=methyl_preprocess

R_SCRIPT=methyl_450_plus_27.R
SH_SCRIPT=methyl_450_plus_27.sh

# Input data
CANCERS=(BLCA BRCA CESC COAD GBM HNSC KIRC KIRP LAML LGG LIHC LUAD LUSC OV PRAD SARC SKCM STAD THCA UCEC)

cd $MY_HOME/$WORKING_DIR

for i in {0..19}; do
	echo $i-${CANCERS[i]}-start qsub
	qsub -v MY_HOME=$MY_HOME,WORKING_DIR=$WORKING_DIR,R_SCRIPT=$R_SCRIPT,CANCER=${CANCERS[i]} $SH_SCRIPT
    echo $i-${CANCERS[i]}-end qsub
done
