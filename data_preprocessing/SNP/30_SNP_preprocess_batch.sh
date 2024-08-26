#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/
WORKING_DIR=$HOME_DIR/scripts/SNP_preprocess/

# Script Description
JOB=30_SNP_preprocess
R_SCRIPT=$WORKING_DIR/${JOB}.R

OE=$WORKING_DIR/OE_$JOB/
if [[ -d "$OE" ]]; then rm -r $OE; fi
mkdir -p $OE

# separate model
DATASETS=(READ-COAD LIHC BLCA KIRC PAAD SARC THCA BRCA ESCA KIRP LUAD PCPG SKCM CESC GBM LUSC PRAD STAD UCEC COAD HNSC LGG OV READ TGCT)
# second run Sep 10
DATASETS=(ESCA LUSC KIRC)
# joint model
# DATASETS=(AVG)

for i in {0..2}; do # todo: chagne when need

    DATASET=${DATASETS[i]}
    INPUT_DIR=$WORKING_DIR/temp/$DATASET/

    OUTPUT_DIR=$HOME_DIR/ref/target_matrix/
    mkdir -p $OUTPUT_DIR

    for j in {0..0}; do

        if [[ $j -eq 0 ]]; then
            EXON_SPECIFIC=$j
            SNP_TARGET_MATRIX=$OUTPUT_DIR/30_SNP_target_matrix_gene_$DATASET.rds
        else
            EXON_SPECIFIC=$j
            SNP_TARGET_MATRIX=$OUTPUT_DIR/30_SNP_target_matrix_exon_$DATASET.rds
        fi

        R_JOB=" module load R/3.3.2;

                Rscript --vanilla $R_SCRIPT
                        --INPUT_DIR $INPUT_DIR
                        --EXON_SPECIFIC $EXON_SPECIFIC
                        --SNP_TARGET_MATRIX $SNP_TARGET_MATRIX"

        echo $R_JOB  qsub -N ${DATASET}-${JOB}-${EXON_SPECIFIC} \
                             -l mem=100G,vmem=100G,walltime=23:59:59 \
                             -o $OE -e $OE

        echo $R_JOB | qsub -N ${DATASET}-${JOB}-${EXON_SPECIFIC} \
                             -l mem=100G,vmem=100G,walltime=23:59:59 \
                             -o $OE -e $OE
    done
done

