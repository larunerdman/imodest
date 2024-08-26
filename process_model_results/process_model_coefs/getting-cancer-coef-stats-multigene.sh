#!/bin/bash

echo "Cancer: " $CANCER 
echo "WORKING_DIR: " $WORKING_DIR

module load R/3.4.4

for REG in {6..6} ; do echo "Running cancer: " $CANCER " and Reg: " $REG ; Rscript --vanilla $WORKING_DIR/get.mean.coefs-function-multigene.R --CANCER $CANCER --REGULATORS $REG ; done




