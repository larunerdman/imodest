#!/bin/bash

#PBS -l vmem=8G,walltime=8:00:00

# loading
module load R/3.2.3

# Dependency

## Main()
#########
echo $R_SCRIPT-----$CANCER start

data_folder=$MY_HOME/$CANCER/methylation_preprocess

cd $data_folder

cp $MY_HOME/$WORKING_DIR/$R_SCRIPT .

Rscript --vanilla ./$R_SCRIPT -A $CANCER

rm ./$R_SCRIPT

echo $R_SCRIPT-----$CANCER finish