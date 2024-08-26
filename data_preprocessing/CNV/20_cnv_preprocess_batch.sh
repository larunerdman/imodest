#!/bin/bash

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/cnv_preprocess/

JOB=20_cnv_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R


CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)


for i in {0..28}; do

    CANCER=${CANCERS[i]}

    # OE
    OE=$WORKING_DIR/OE_$JOB/$CANCER/
    if [[ -d "$OE" ]]; then rm -r "$OE"; fi
    mkdir -p "$OE"

    # Input specification
    CNV_MATRIX_RAW=$HOME_DIR/cnv_preprocess_output/raw/${CANCER}_cnv_matrix_gene.rds

    #CLINICAL_DATA=$HOME_DIR/clinical_preprocess_output/${CANCER}_clinical_data.rds
    CLINICAL_DATA=/hpf/largeprojects/agoldenb/mingjie_new/clinical_preprocess_output/${CANCER}_clinical_data.rds

    # Output specification
    #OUTPUT_DIR=$HOME_DIR/cnv_preprocess_output/corrected/$CANCER/
    OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/cnv_preprocess_output/corrected/$CANCER/
    TEMP_DIR=$OUTPUT_DIR/correct_data_analysis/
    if [[ -d "$OUTPUT_DIR" ]]; then rm -r "$OUTPUT_DIR"; fi
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

    CNV_MATRIX_CORRECTED=$OUTPUT_DIR/${CANCER}_cnv_matrix_gene_corrected.rds

    # Other
    HELPER=$HOME_DIR/helper/correct_data.R

    NUM_CORES=2

    R_JOB="module load R/3.3.2;

            Rscript --vanilla $R_SCRIPT
                    --CANCER $CANCER
                    --CNV_MATRIX_RAW $CNV_MATRIX_RAW
                    --CLINICAL_DATA $CLINICAL_DATA
                    --CNV_MATRIX_CORRECTED $CNV_MATRIX_CORRECTED
                    --TEMP_DIR $TEMP_DIR
                    --HELPER $HELPER
                    --NUM_CORES $NUM_CORES"

    echo $R_JOB   qsub -N $CANCER-$JOB \
                       -l mem=30G,vmem=30G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                       -p 0 \
                       -o $OE -e $OE

    echo $R_JOB | qsub -N $CANCER-$JOB \
                       -l mem=30G,vmem=30G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                       -p 0 \
                       -o $OE -e $OE
done
