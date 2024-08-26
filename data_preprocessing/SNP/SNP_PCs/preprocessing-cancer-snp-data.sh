#!/bin/bash

#PBS -l vmem=30g,walltime=23:00:00


(>&2 echo WORKING_DIR: $WORKING_DIR
>&2 echo CANCER: $CANCER
>&2 echo DATA_DIR: $DATA_DIR)


cd $DATA_DIR

awk '{if(! a[$2]){print ; a[$2]++}}' $CANCER.tped | grep -v '-' > $CANCER-no-.tped

module load plink
module load python/2.7.9

python $WORKING_DIR/preprocessing-genome-wide-snp-data-id_sub.py -c $CANCER
