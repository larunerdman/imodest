#!/bin/bash

umask 0007

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/RNA_preprocess/

#### Job Description
JOB=20_RNA_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R


CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)
#CANCERS=(OV)
COME_FROM_TAR=0

PEER_HIDDENS_SAMPLES_RATIOS=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95)

for i in {0..28}; do

    CANCER=${CANCERS[i]}

    # Make OE
    OE=$WORKING_DIR/OE_$JOB/peer/$CANCER/
    if [[ -d "$OE" ]]; then rm -r "$OE"; fi
    mkdir -p "$OE"

    # Input specification
    CLINICAL_DATA_RAW=$HOME_DIR/clinical_preprocess_input/${CANCER}.clin.merged.txt
    RNA_INPUT_FOLDER=$HOME_DIR/RNA_preprocess_input/${CANCER}/gene/

    #CLINICAL_DATA=$HOME_DIR/clinical_preprocess_output/${CANCER}_clinical_data.rds
    CLINICAL_DATA=/hpf/largeprojects/agoldenb/mingjie_new/clinical_preprocess_output/${CANCER}_clinical_data.rds
    TF_TARGET_MATRIX=$HOME_DIR/ref/target_matrix/10_TF_target_matrix_gene.rds

    # Output specification
    #OUTPUT_DIR=$HOME_DIR/RNA_preprocess_output/peer/$CANCER/
    OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/RNA_preprocess_output/peer/$CANCER/
    TEMP_DIR=$OUTPUT_DIR/correct_data_analysis/
    # if [[ -d "$OUTPUT_DIR" ]]; then rm -r "$OUTPUT_DIR"; fi
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

    # Other
    HELPER=$HOME_DIR/helper/correct_data.R

    NUM_CORES=4

    for j in {1..19}; do

        PEER_HIDDENS_SAMPLES_RATIO=${PEER_HIDDENS_SAMPLES_RATIOS[j]}

        PC_GENE_EXPRESSION_CORRECTED=$OUTPUT_DIR/${CANCER}_pc_gene_expression_corrected_$j.rds
        LNCRNA_GENE_EXPRESSION_CORRECTED=$OUTPUT_DIR/${CANCER}_lncRNA_gene_expression_corrected_$j.rds
        TF_GENE_EXPRESSION_CORRECTED=$OUTPUT_DIR/${CANCER}_TF_gene_expression_corrected_$j.rds

        #if [[ ! -f $PC_GENE_EXPRESSION_CORRECTED || ! -f $LNCRNA_GENE_EXPRESSION_CORRECTED || ! -f $TF_GENE_EXPRESSION_CORRECTED ]]; then
            R_JOB=" module load R/3.3.2;

                    Rscript --vanilla $R_SCRIPT
                            --CANCER $CANCER
                            --CLINICAL_DATA_RAW $CLINICAL_DATA_RAW
                            --RNA_INPUT_FOLDER $RNA_INPUT_FOLDER
                            --CLINICAL_DATA $CLINICAL_DATA
                            --TF_TARGET_MATRIX $TF_TARGET_MATRIX
                            --TEMP_DIR $TEMP_DIR
                            --PEER_HIDDENS_SAMPLES_RATIO $PEER_HIDDENS_SAMPLES_RATIO
                            --PC_GENE_EXPRESSION_CORRECTED $PC_GENE_EXPRESSION_CORRECTED
                            --LNCRNA_GENE_EXPRESSION_CORRECTED $LNCRNA_GENE_EXPRESSION_CORRECTED
                            --TF_GENE_EXPRESSION_CORRECTED $TF_GENE_EXPRESSION_CORRECTED
                            --HELPER $HELPER
                            --NUM_CORES $NUM_CORES
                            --COME_FROM_TAR $COME_FROM_TAR"

            echo $R_JOB   qsub -N $CANCER-$JOB \
                               -l mem=60G,vmem=60G,walltime=100:00:00,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                               -p 0 \
                               -o $OE -e $OE

            echo $R_JOB | qsub -N $CANCER-$JOB \
                               -l mem=60G,vmem=60G,walltime=100:00:00,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                               -p 0 \
                               -o $OE -e $OE
       #fi
   done
done






