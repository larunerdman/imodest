#!/bin/bash

#PBS -l vmem=30g,walltime=23:00:00


(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo PREFIX: $PREFIX
>&2 echo DATA_DIR: $DATA_DIR)


cd $DATA_DIR

module load plink
module load python/2.7.9

python $WORKING_DIR/preprocessing-genome-wide-snp-data-nonTCGA.py -p $CANCER
