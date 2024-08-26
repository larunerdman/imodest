#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/methylation_preprocess/

JOB=50_methylation_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R

OE=$WORKING_DIR/OE_$JOB/
mkdir -p $OE
if [ "$(ls -A $OE)" ]; then rm $OE/*; fi

for j in {1..1}; do

    INPUT_DIR=$WORKING_DIR/temp/

    EXON_SPECIFIC=$j

    if [[ $EXON_SPECIFIC -eq 0 ]]; then
        METHYLATION_TARGET_MATRIX=$HOME_DIR/ref/target_matrix/50_methylation_target_matrix_gene.rds
    else
        METHYLATION_TARGET_MATRIX=$HOME_DIR/ref/target_matrix/50_methylation_target_matrix_exon.rds
    fi

    R_JOB=" module load R/3.3.2;

            Rscript --vanilla $R_SCRIPT
                    --INPUT_DIR $INPUT_DIR
                    --EXON_SPECIFIC $EXON_SPECIFIC
                    --METHYLATION_TARGET_MATRIX $METHYLATION_TARGET_MATRIX"

    echo $R_JOB | qsub -N $JOB-$EXON_SPECIFIC \
                       -l vmem=60G,walltime=23:59:59 \
                       -o $OE -e $OE
done


