#!/bin/bash

#PBS -l vmem=2G,walltime=23:59:59
umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/projects/MoR/
WORKING_DIR=$HOME_DIR/scripts/SNP_preprocess/

JOB=40_SNP_preprocess

R_SCRIPT=$WORKING_DIR/$JOB.R


# separate model
DATASETS=(READ-COAD LIHC BLCA KIRC PAAD SARC THCA BRCA ESCA KIRP LUAD PCPG SKCM CESC GBM LUSC PRAD STAD UCEC COAD HNSC LGG OV READ TGCT)
# second run Sep 10
DATASETS=(ESCA LUSC KIRC)
# joint model
# DATASETS=(AVG)

for i in {0..2}; do # todo: chagne when need

    DATASET=${DATASETS[i]}

    OE=$WORKING_DIR/OE_$JOB/$DATASET/
    if [[ -d "$OE" ]]; then rm -r $OE; fi
    mkdir -p $OE

    # Input folder (imputed SNP)
    # SNP_INPUT_DIR=$HOME_DIR/SNP_preprocess_input/$DATASET/genotype-calls-july-2016/imputation/
    SNP_INPUT_DIR=$HOME_DIR/data/individual_cohort/SNP_preprocess_input/$DATASET/genotype-calls-july-2016/imputation/

    for j in {1..22}; do

        # Input
        SAMPLE_FILE=$SNP_INPUT_DIR/${DATASET}.sample

        cat $SNP_INPUT_DIR/$DATASET-cancer-postqc-postqc2-revised-chr1.sample | grep "TCGA" | cut -c 1-16 > $SAMPLE_FILE
        BIM_FILE=$SNP_INPUT_DIR/$DATASET-cancer-postqc-postqc2-postqc-somup-chr$j.bim
        DOSAGE_FILE=$SNP_INPUT_DIR/$DATASET-cancer-postqc-postqc2-postqc-somup-chr$j.dosage

        # cat $SNP_INPUT_DIR/all-TCGA-mj-sub-postqc2-revised-chr1.sample | grep "TCGA" | cut -c 1-16 > $SAMPLE_FILE
        # BIM_FILE=$SNP_INPUT_DIR/all-TCGA-mj-sub-postqc2-postqc-chr$j.bim
        # DOSAGE_FILE=$SNP_INPUT_DIR/all-TCGA-mj-sub-postqc2-postqc-chr$j.dosage

        # Output
        SNP_OUTPUT_DIR=$HOME_DIR/SNP_preprocess_output/$DATASET/
        if [[ -d "$SNP_OUTPUT_DIR" ]]; then rm -r $SNP_OUTPUT_DIR; fi
        mkdir -p $SNP_OUTPUT_DIR

        SNP_MATRIX=$SNP_OUTPUT_DIR/${DATASET}_SNP_matrix_chr$j.rds

        R_JOB=" module load R/3.3.2;

                Rscript --vanilla $R_SCRIPT
                        --BIM_FILE $BIM_FILE
                        --DOSAGE_FILE $DOSAGE_FILE
                        --SAMPLE_FILE $SAMPLE_FILE
                        --SNP_MATRIX $SNP_MATRIX"

        echo $R_JOB    qsub -N $DATASET-${JOB}-$j \
                             -l mem=60G,vmem=60G,walltime=23:59:59 \
                             -o $OE -e $OE
        echo $R_JOB |   qsub -N $DATASET-${JOB}-$j \
                             -l mem=60G,vmem=60G,walltime=23:59:59 \
                             -o $OE -e $OE
    done
done

