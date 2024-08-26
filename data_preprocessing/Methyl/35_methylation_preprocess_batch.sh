#!/bin/bash

# ------------
# Do standardization, covariate adjustment or peer correction on the raw
#   methylation data
# ------------

HOME_DIR=/hpf/largeprojects/agoldenb/mingjie/
WORKING_DIR=$HOME_DIR/methylation_preprocess/

JOB=35_methylation_preprocess
R_SCRIPT=$WORKING_DIR/$JOB.R


# Input data
CANCERS=(READ-COAD LIHC BLCA COAD BRCA KIRC LGG HNSC LUAD PRAD SKCM UCEC CESC THCA KIRP LUSC PAAD TGCT PCPG GBM THYM SARC READ LAML STAD ESCA OV STAD ESCA)
PEER_HIDDENS_SAMPLES_RATIOS=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95)



for i in {0..28}; do


    CANCER=${CANCERS[i]}

    OE=$WORKING_DIR/OE_$JOB/peer/$CANCER/
    #if [[ -d "$OE" ]]; then rm -r "$OE"; fi
    mkdir -p "$OE"

    # Input specification
    METHYLATION_MATRIX_RAW=$HOME_DIR/methylation_preprocess_output/raw/${CANCER}_methylation_matrix.rds

    #CLINICAL_DATA=$HOME_DIR/clinical_preprocess_output/${CANCER}_clinical_data.rds
    CLINICAL_DATA=/hpf/largeprojects/agoldenb/mingjie_new/clinical_preprocess_output/${CANCER}_clinical_data.rds
    #OUTPUT_DIR=$HOME_DIR/methylation_preprocess_output/peer/$CANCER/
    OUTPUT_DIR=/hpf/largeprojects/agoldenb/mingjie_new/methylation_preprocess_output/peer/$CANCER/
    TEMP_DIR=$OUTPUT_DIR/correct_data_analysis/
    # if [[ -d "$OUTPUT_DIR" ]]; then rm -r "$OUTPUT_DIR"; fi
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

    # Other
    HELPER=$HOME_DIR/helper/correct_data.R

    NUM_CORES=4

    for j in {0..19}; do

        for d in {1..2}; do
            DATA_SET=$d
            PEER_HIDDENS_SAMPLES_RATIO=${PEER_HIDDENS_SAMPLES_RATIOS[j]}
            METHYLATION_MATRIX_CORRECTED=$TEMP_DIR/${CANCER}_methylation_matrix_corrected_${j}_${d}.rds

            #if [[ ! -f $METHYLATION_MATRIX_CORRECTED ]]; then
                R_JOB=" module load R/3.3.2;

                        Rscript --vanilla $R_SCRIPT
                                --CANCER $CANCER
                                --METHYLATION_MATRIX_RAW $METHYLATION_MATRIX_RAW
                                --CLINICAL_DATA $CLINICAL_DATA
                                --METHYLATION_MATRIX_CORRECTED $METHYLATION_MATRIX_CORRECTED
                                --TEMP_DIR $TEMP_DIR
                                --PEER_HIDDENS_SAMPLES_RATIO $PEER_HIDDENS_SAMPLES_RATIO
                                --HELPER $HELPER
                                --DATA_SET $DATA_SET
                                --NUM_CORES $NUM_CORES"

                echo $R_JOB   qsub -N $CANCER-$JOB \
                                   -l mem=60G,vmem=60G,walltime=100:00:00,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                                   -p 100 \
                                   -o $OE -e $OE

                echo $R_JOB | qsub -N $CANCER-$JOB \
                                   -l mem=60G,vmem=60G,walltime=100:00:00,nodes=1:ppn=$NUM_CORES,gres=localhd:0 \
                                   -p 100 \
                                   -o $OE -e $OE
            #fi
        done

    done

done
