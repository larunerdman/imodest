#!/bin/bash

#PBS -l vmem=1G,walltime=1:00:00
HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/miRNA_preprocess/

JOB=20_miRNA_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R


HELPER=$HOME_DIR/helper/correct_data.R

# Input data
CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)
PEER_HIDDENS_SAMPLES_RATIOS=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95)

for i in {0..28}; do

    CANCER=${CANCERS[i]}

    # Make OE
    OE=$WORKING_DIR/OE_$JOB/$CANCER/
    if [[ -d "$OE" ]]; then rm -r "$OE"; fi
    mkdir -p "$OE"


    # Input specification
    MIRNA_EXPRESSION_DATA=$HOME_DIR/miRNA_preprocess_input/raw/${CANCER}.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data.data.txt
    #CLINICAL_DATA=$HOME_DIR/clinical_preprocess_output/${CANCER}_clinical_data.rds
    CLINICAL_DATA=/hpf/largeprojects/agoldenb/mingjie_new/clinical_preprocess_output/${CANCER}_clinical_data.rds
    NUM_CORES=4


    # Output specification
    OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/miRNA_preprocess_output/peer/$CANCER/
    TEMP_DIR=$OUTPUT_DIR/correct_data_analysis/
    #if [[ -d "$OUTPUT_DIR" ]]; then rm -r "$OUTPUT_DIR"; fi
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

    for j in {1..19}; do
        PEER_HIDDENS_SAMPLES_RATIO=${PEER_HIDDENS_SAMPLES_RATIOS[j]}
        MIRNA_EXPRESSION_MATRIX_CORRECTED=$OUTPUT_DIR/${CANCER}_miRNA_expression_matrix_corrected_$j.rds
        R_JOB=" module load R/3.3.2; \
        	    Rscript --vanilla $R_SCRIPT
                	    --CANCER $CANCER
                        --MIRNA_EXPRESSION_DATA $MIRNA_EXPRESSION_DATA
                        --CLINICAL_DATA $CLINICAL_DATA
                        --TEMP_DIR $TEMP_DIR
                        --MIRNA_EXPRESSION_MATRIX_CORRECTED $MIRNA_EXPRESSION_MATRIX_CORRECTED
                        --PEER_HIDDENS_SAMPLES_RATIO $PEER_HIDDENS_SAMPLES_RATIO
                        --HELPER $HELPER
                        --NUM_CORES $NUM_CORES"

    	echo $R_JOB   qsub -N $CANCER-$JOB \
                           -l mem=60G,vmem=60G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                           -p 0 \
                           -o $OE -e $OE

        echo $R_JOB | qsub -N $CANCER-$JOB \
                           -l mem=60G,vmem=60G,walltime=23:59:59,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                           -p 0 \
                           -o $OE -e $OE
    done
done
