#!/bin/bash

#PBS -l vmem=10g,walltime=15:00:00

echo PARFILE: $PARFILE
echo LOGFILE: $LOGFILE
echo DATA_DIR: $DATA_DIR

cd $DATA_DIR

module load eigensoft/6.0

smartpca -p $PARFILE > $LOGFILE
